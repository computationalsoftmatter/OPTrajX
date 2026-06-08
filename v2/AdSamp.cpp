/*License:
--------
This software is distributed under the GNU General Public License (GPL),
either version 3 of the License or (at your option) any later version.
You are free to use, modify, and distribute this software under the terms
of the GPL. For more details, see the LICENSE file included with this project
or visit https://www.gnu.org/licenses/gpl-3.0.html.

Contributors:
-------------
- Manuel Lopez
- Betul Uralcan
- Amir Haji-akbari

Contact:
--------
For questions, please contact:
- Betul Uralcan
  - Email: betul.uralcan@bogazici.edu.tr
  - Email: betul.uralcan@yale.edu

*/

#include <filesystem>
#include <iostream>
#include <chrono>
#include <thread>
#ifdef USE_LAMMPS
#include "lammps/lammps.h"
#include "lammps/input.h"
#include "lammps/library.h"
#include "lammps/atom.h"
#endif
#include "mpi.h"
#include <sstream>
#include <vector>
#include <random>
#include <fstream>
#include <algorithm>
#include "orderparameterpicker.hpp"
#include <cstdlib>
#include <cstdio>
#include "OrderParameter.hpp"
#include <regex>
#include <string>
#include <cctype>
#include <ctime>
#include <unordered_map> 
#include <csignal>
#include <atomic>

volatile sig_atomic_t stop_requested = 0;

void signal_handler(int signum) {
    stop_requested = signum; //for lammps interruptions
}

//alright double triple check this parallelization stuff
struct Context { //okay the idea here is to have a serial context for gromacs and a parallel one for LAMMPS
    virtual ~Context()=default;
    virtual int rank() const = 0;
    virtual int size() const = 0;
    virtual void barrier() const = 0;
    virtual void bcast_int(int&x , int root) const=0;//have to make these real for parallel and fake for series
    virtual void bcast_longlong(long long& x, int root ) const =0;
    virtual void bcast_double(double&x, int root) const = 0;
    virtual void bcast_string(std::string& s, int root) const  = 0;
    virtual MPI_Comm get_comm() const = 0; 
};

struct ParallelContext: Context {
    MPI_Comm comm; //wowza so the explanation for this stuff and it's need is a little long winded
    int r; //long story short lammps needs us to split into paralel cores but if we do that in gromacs we;d only be able to run the actual sim on one core or do it n times
    int n;  //the solution is to split if we're on lammps but not othewrise, so to do this we mimic the way we insert an order parameter object with a virtual class
    explicit ParallelContext(MPI_Comm c=MPI_COMM_WORLD) : comm(c) { //so it's like a virtual MPI handler which is conceptually pretty cool
        MPI_Comm_rank(comm, &r);
        MPI_Comm_size(comm, &n);
    }
    int rank() const override { 
        return r; 
    }
    int size() const override { 
        return n; 
    }
    void bcast_longlong(long long& x, int root) const override {
        MPI_Bcast(&x, 1, MPI_LONG_LONG, root, comm);
    }
    void barrier() const override { MPI_Barrier(comm); }
    void bcast_int(int&x , int root) const override{
        MPI_Bcast(&x, 1, MPI_INT, root, comm);
    }
    void bcast_double(double&x, int root) const override{
        MPI_Bcast(&x, 1, MPI_DOUBLE, root, comm);
    }
    void bcast_string(std::string &s, int root) const override {
        int myrank = -1;
        MPI_Comm_rank(comm, &myrank);

        int len = (myrank == root) ? static_cast<int>(s.size()) : 0;
        MPI_Bcast(&len, 1, MPI_INT, root, comm);

        if (myrank != root) s.resize(len);

        if (len > 0) {
            MPI_Bcast(&s[0], len, MPI_CHAR, root, comm);
        }
    }
    MPI_Comm get_comm() const override { return comm; }
};

struct SerialContext : Context{//okay this ones way simpler, we always return rank 0 and size 1, no barrier, no broadcast
    int rank() const override{ 
        return 0;
    }
    int size() const override{
        return 1;
    }
    void bcast_longlong(long long& x, int root) const override { (void)x; (void)root; }
    void barrier() const override{

    }
    void bcast_int(int& x, int root) const override { (void)x; (void)root; }
    void bcast_double(double& x, int root) const override { (void)x; (void)root; }
    void bcast_string(std::string& s, int root)  const override { 
        (void)s; (void)root; 
    } 
    MPI_Comm get_comm() const override { return MPI_COMM_SELF; }
};


static std::string lastdataline(const std::string& colvar) { //faster way to get the last line of a colvar file, used as a lammps speedboost
    std::ifstream in(colvar, std::ios::binary | std::ios::ate);
    if (!in) return "";
    std::streamoff pos = static_cast<std::streamoff>(in.tellg());
    if (pos == 0) return "";
    std::string line;
    while (pos > 0) {
        in.seekg(--pos);
        char c = static_cast<char>(in.peek());
        if (c != '\n' && c != '\r' && c != ' ') break;
    }
    while (pos > 0) {
        in.seekg(--pos);
        char c = static_cast<char>(in.peek());
        if (c == '\n' || c == '\r') { ++pos; break; }
    }
    in.seekg(pos);
    std::getline(in, line);
    if (line.empty() || line[0] == '#') return "";
    return line;
}


std::vector<std::string> readlines(const std::string& path) {
    std::ifstream in(path);
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(in, line)){
        lines.push_back(line);
    }
    return lines;
}

void DestroyLammps(LAMMPS_NS::LAMMPS*& lmp)
{
    if (lmp) {
        delete lmp;
        lmp = nullptr;  
    }
}
std::vector<std::string> tokenize(std::string line){
    std::istringstream iss(line);
    std::vector<std::string> tokens;
    std::string token;
    while (iss>>token){
        tokens.push_back(token);
    }
    return tokens;
}
double ReadTemp(std::string path){
    std::ifstream mdp(path);
    std::string line;
    while (std::getline(mdp, line)){
        std::vector<std::string> tokens = tokenize(line);
        for (int i=0; i< tokens.size(); i++ ){
            if (tokens[i] =="npt" || tokens[i] =="nvt"){ //temperature line must have one of these
                for (int j=i+1; j<tokens.size(); j++){
                    if (tokens[j] =="temp" && j< tokens.size()-1){
                        return std::stod(tokens[j+1]) ;//
                    }
                }

            }
        }
    }
    return 1.0; //shouldn't ever need this but to compile
}
std::pair<double,double> lastorderparameter(std::string colvar){ //so truthfully this need not be its own function, however, there's a good chance that LAMMPS and plumed print op in very different ways, this will help provide some uniformity
    std::string lastline = lastdataline(colvar);
    std::stringstream stream(lastline);
    double time, value;
    stream>>time>>value;
    return {time, value};
}

int computecriteria(long long paso, long long nsteps, double orderparam,  int direction, double leftbound, double rightbound, int runtype){ //this function implements the 'handedness' of basin termination criteria
    if (paso>=nsteps){
        return 1;// okay so we should never be greater than but better safe than sorry
    }
    if (runtype ==1){ //right handed FFS
        if (orderparam*direction >= direction*rightbound){
            return 1;
        }else{
            return 0;
        }
    }
    if (runtype ==2){ //left handed FFS
        if (orderparam*direction<= direction*leftbound){
            return 1;
        }else{
            return 0;
        }
    }if (runtype ==3){//standard FFS for probabilistic step
        if (orderparam*direction<= direction*leftbound || orderparam*direction>=direction*rightbound){
            return 1;
        }else{
            return 0;
        }
    }// note that in this workflow I don't have a condition for runtype 4 other than nsteps
    return 0;
}
void Colvar_Writer(std::string nombre, double time, double OP, double deltaT){
    std::ofstream out(nombre, std::ios::app);
    out<< time*deltaT<< " "<< OP<< " "<< std::endl;
}

#ifdef USE_LAMMPS //all this uselammps stuff is so that a person can compile and run this for gromacs without also compiling LAMMPS
void LAMMPSVelocities(LAMMPS_NS::LAMMPS* lmp, double Temp){
    //note first go already sets random seed
    std::normal_distribution<double> gauss(0.0,1.0); //make them doubles (0.0)
    std::mt19937 rng(std::random_device{}()); //create random number generator
    double** v=lmp->atom->v;
    double* m=lmp->atom->mass;
    int* type=lmp->atom->type;
    int nlocal=lmp->atom->nlocal;
    double kT= Temp; //boltzmann factor is one for standard lammps units
    for (int i =0; i<nlocal; i++){
        double scale = sqrt( kT/m[type[i]]);
        v[i][0]=gauss(rng)*scale;
        v[i][1]=gauss(rng)*scale;
        v[i][2]=gauss(rng)*scale;
    }
}

LAMMPS_NS::LAMMPS* InitializeLammps(Context* ctx,std::string newfolder, std::string propername, std::string structurefile, int feathered ){//little plumed wordplay over here
    //I read  dummy values are fine, come back to check this later
    bool isCPT = structurefile.size()>=4 && structurefile.substr(structurefile.size()-4 )==".cpt"; //okay so if we're in the lammps workflow and the second run uses the checkpoint of a prior run then we need to do a diferent command as its binary
    std::string prelim = isCPT ? "read_restart " : "include "; 
    int argc = 0;
    char **argv = nullptr;
    LAMMPS_NS::LAMMPS* lmp = new LAMMPS_NS::LAMMPS(argc, argv, ctx->get_comm());
    std::string cmd = prelim + structurefile; //structure file command
    lmp->input->one(cmd.c_str());
    if (isCPT){
        lmp->input->one("reset_timestep 0"); //hopefully this works, lammps uses the cpt's time naturally so it doesn't go back to zero
    }
    cmd = "include " + newfolder + "/" + propername + ".mdp";
    lmp->input->one(cmd.c_str());
    double Temp = ReadTemp(newfolder + "/" + propername + ".mdp");
    LAMMPSVelocities(lmp,Temp);
    if (feathered ==1){
        cmd = "fix cv all plumed plumedfile "  + newfolder +"/" + propername +".dat outfile " + newfolder + "/" + propername + ".log";
        lmp->input->one(cmd.c_str());
    }
    return lmp;
}
LAMMPS_NS::LAMMPS* RestartLammps(Context* ctx,std::string newfolder, std::string propername, int feathered){
    int argc=0;
    char **argv =nullptr;
    LAMMPS_NS::LAMMPS* lmp = new LAMMPS_NS::LAMMPS(argc, argv, ctx->get_comm());
    std::string cmd="read_restart " + newfolder +"/" + propername+".cpt";
    lmp->input->one(cmd.c_str());
    lmp->input->one("reset_timestep 0"); 
    cmd ="include " + newfolder + "/" + propername + ".mdp";
    lmp->input->one(cmd.c_str());
    if (feathered ==1){
        cmd = "fix cv all plumed plumedfile " + newfolder + "/" + propername +".dat outfile " + newfolder + "/" + propername + ".log";
        lmp->input->one(cmd.c_str());
    }
    return lmp;
}
void writeCPT(LAMMPS_NS::LAMMPS* lmp, std::string newfolder, std::string propername){
    std::string restartFile = newfolder + "/" + propername + ".cpt"; 
    if (std::filesystem::exists(restartFile)){
        std::filesystem::remove(restartFile);
    }
    std::string cmd="write_restart " + restartFile;
    lmp->input->one(cmd.c_str()); 
}
void lammpsrunner(Context* ctx,LAMMPS_NS::LAMMPS* lmp, OrderParameter* op, int usePLMD, long long nsteps, long long runtime, int direction, int runtype, double leftbound, double rightbound,std::string newfolder, std::string propername, int strides, double deltaT){
    const int rank = ctx->rank();
    stop_requested = 0;
    std::signal(SIGUSR1, signal_handler);
    std::signal(SIGTERM, signal_handler);
    std::signal(SIGUSR2, signal_handler);
    int exit = 0;

    if (usePLMD ==1){
        long long paso=0;
        double time;
        double c;
        std::string RColvarPath = newfolder + "/RC" + propername;
        std::string ColvarPath = newfolder + "/C" + propername;
        long long remaining = nsteps - runtime;
        long long largest = strides; //break it down into stride size steps for cpt safety
        bool first = true;
        while (remaining>0){ //this architecture is so that a user can use nsteps larger than what lammps accepts as max steps
            long long now =std::min(remaining, largest);
            std::string correr = "run " + std::to_string(now) + (first ? " pre yes post no" : " pre no post no");
            lmp->input->one(correr.c_str());
            if (stop_requested){
                writeCPT(lmp, newfolder, propername);
                return;
            }
            first = false;
            remaining -=now;

            if (rank ==0 && runtype !=4){
                std::string which;
                int use = 0;
                if (std::filesystem::exists(RColvarPath)){
                    which = RColvarPath;
                    use = 1;
                }else if(std::filesystem::exists(ColvarPath)){
                    which = ColvarPath;
                    use =1;
                }
                if (use ==1){
                    auto result = lastorderparameter(which);
                    time = result.first; //no need to convert to nsteps here
                    paso = time/deltaT;
                    c= result.second;
                    std::cout<<c<<std::endl;
                    int done = computecriteria(paso,nsteps, c, direction, leftbound, rightbound, runtype);
                    if (done){
                        exit =1;
                    }
                }
            }
            ctx->barrier();
            ctx->bcast_int(exit, 0);
            if (exit==1){
                return;
            }
        }
    }else{
        std::cout<<"runtype: "<< runtype<<std::endl;
        std::string prelim = "/C";
        if (std::filesystem::exists(newfolder + "/C" + propername)){
            prelim = "/RC";
        }
        
        long long paso = runtime;
        int micondicion = 0;
        while (micondicion ==0){
            long long now = std::min((long long)strides, nsteps - paso);
            std::string cmd = "run " +std::to_string(now)+ " pre no post no";
            lmp->input->one(cmd.c_str());
            if (stop_requested){
                writeCPT(lmp, newfolder, propername);
                return;
            }
            paso += now; 
            double orderparam = op->update(lmp);
            if (rank ==0){
                Colvar_Writer(newfolder + prelim +propername, paso, orderparam ,deltaT);
            }
            ctx->barrier();
            std::cout<<"Order Parameter "<< op<<std::endl;
            micondicion = computecriteria(paso, nsteps, orderparam, direction, leftbound, rightbound, runtype);
            std::cout<<micondicion<<std::endl;
        }
    }
}

#endif






int getenv_int(const char* name, int default_value) {
    const char* val = std::getenv(name);
    if (!val) return default_value;
    return std::stoi(val);
}
static int slurm_total_tasks() {
    return getenv_int("SLURM_NTASKS", 1);
}
void InitializeGMX(std::string MDPFILE,  std::string GROFILE, std::string TOPFILE,std::string NDXFILE, std::string TPRFILE, std::string CPTFILE,Context* ctx ){ //make TPR file
    int rank = ctx->rank();
    if (rank!=0){
        return;
    }
    int numcores = slurm_total_tasks();
    //std::string cmd = "srun --overlap --overcommit --cpu-bind=none --mem-bind=none -n " + std::to_string(numcores) + " gmx_mpi grompp -f " + MDPFILE + " -c " + GROFILE + " -r " + GROFILE + " -p " + TOPFILE + " -n " + NDXFILE + " -o " +TPRFILE + " -maxwarn 2";
    std::string cmd = "srun --overlap --overcommit --cpu-bind=none --mem-bind=none --mem=0 -n 1 gmx_mpi grompp -f " + MDPFILE + " -c " + GROFILE + " -r " + GROFILE + " -p " + TOPFILE + " -n " + NDXFILE + " -o " +TPRFILE + " -maxwarn 2";
    if (!CPTFILE.empty()){
        cmd += " -t " + CPTFILE;
    }
    std::system(cmd.c_str());
}

void RUNGMX(std::string nom, int tipo, Context* ctx ){
    int rank = ctx->rank();
    if (rank!=0){
        return;
    }
    int numcores = slurm_total_tasks();
    std::cout<<"cores"<<numcores<<std::endl;
    //std::string cmd = "srun --overlap --overcommit --cpu-bind=none --mem-bind=none --mem=0 -n " +std::to_string(numcores) + " gmx_mpi mdrun -v -deffnm "+ nom + " -plumed " + nom + ".dat -c " + nom +".gro -append -cpt 5 -ntomp 1";
    //std::string cmd = "srun --overlap --overcommit --cpu-bind=cores --mem-bind=local --exact --mem=0 -n " +std::to_string(numcores) + " gmx_mpi mdrun -v -deffnm "+ nom + " -plumed " + nom + ".dat -c " + nom +".gro -append -cpt 5 -ntomp 1";
    //std::string cmd = "srun --overlap --cpu-bind=cores --mem-bind=local --exact --mem=0 -n " +std::to_string(numcores) + " gmx_mpi mdrun -v -deffnm "+ nom + " -plumed " + nom + ".dat -c " + nom +".gro -append -cpt 5 -ntomp 1";
    std::string cmd = "srun --cpu-bind=none -n " +std::to_string(numcores) + " gmx_mpi mdrun -v -deffnm "+ nom + " -plumed " + nom + ".dat -c " + nom +".gro -append -cpt 5 -ntomp 1";
    if (tipo ==1){
        cmd += " -cpi " + nom + ".cpt";
    }
    int ret = std::system(cmd.c_str()); //this will run gromacs
    std::cout << "[InitializeGMX] grompp exited with code: " << ret << std::endl;
    std::cout.flush();
}
std::string trim(const std::string& s){//whitespace remover
    size_t start = s.find_first_not_of(" \t\r\n");
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end-start+1);
}
std::vector<std::string> pathlinereader(std::string Line){
    std::stringstream ss(Line);
    std::string entry;
    std::vector<std::string> output;
    while (std::getline(ss, entry, ',')){// comma delimitor is the way I was doing it in gromacs
        output.push_back(trim(entry));
    }
    return output;
}

std::vector<std::string> listdirStart(std::string path, std::string start){
    std::vector<std::string> files;
    for (const auto &entry : std::filesystem::directory_iterator(path)){
        std::string name = entry.path().filename().string();
        if (name.rfind(start, 0)== 0){
            files.push_back(entry.path().filename().string()); //.filename() to remove ath
        }
    }
    return files;
}
std::vector<std::string> listdirEnding(std::string path, std::string ending){
    std::vector<std::string> files;
    for (const auto &entry : std::filesystem::directory_iterator(path)){
        if (entry.path().extension() == ending){
            files.push_back( entry.path().filename().string());
        }
    }
    std::sort(files.begin(), files.end()); //added this in a little later, this gurrantees the files are returned the same way every time (iteration coordination)
    return files;
}
int crossings_checker(std::string colvarfile, double leftbound, double rightbound, int direction){
    //std::vector<std::string> colvarcontent = readlines(colvar) man, I wish but would iterate through later so that's inefficient
    std::ifstream in(colvarfile);
    std::string line;
    double time, op;
    int crossings=0;
    int returned = 1;
    while (std::getline(in, line)){
        if (line.empty() || line[0] == '#') continue; //skip header
        std::istringstream iss(line); //extracting numbers from files in c++ is a little weird
        iss>>time>>op;
        if (direction*op>= direction*rightbound && returned ==1){
            crossings++;
            returned = 0;
        }
        if (direction*op<=direction*leftbound){
            returned =1;
        }
    }
    return crossings;

}
std::string changenumber( std::string line, double value){
    std::size_t spot = line.find("=");
    std::ostringstream output;
    output << line.substr(0, spot + 1) << value; 
    return output.str();
}
void RCfilecombiner(std::string newfolder, std::string RCfile, double PLUMEDSTEPSIZE, int reducesize){
    std::string Cfile = RCfile.substr(1);
    if (!std::filesystem::exists(newfolder + "/"+Cfile)){
        std::filesystem::rename(newfolder + "/" + RCfile,newfolder + "/" + Cfile); //if some glitch happened and we got rid of the C file, kinda a failsafe, I saw it happen once with diskquota
    }else{
        std::vector<std::string> colvarcontent = readlines(newfolder +"/" + Cfile);
        std::vector<std::string> rcolvarcontent =readlines(newfolder +"/" + RCfile);
        if (rcolvarcontent.size() < reducesize+1){
            std::filesystem::remove(newfolder + "/" + RCfile);
            return; //in case rcolvar is empty or otherwise not useful
        }
        //first we gotta get what the las ttime was from the colvar, then we start merging but add this value to the rc values, let's make this a global merger function
        double LastTime;
        std::string lastline= colvarcontent.back();
        std::stringstream stream(lastline);
        stream>>LastTime;
        // let's now find out the first time of the RC
        double firstTime;
        std::string RFirst =rcolvarcontent[reducesize];
        std::stringstream stream1(RFirst);
        stream1>>firstTime;
        //let's start building the merged list
        std::vector<std::string> mergedlist;
        mergedlist.reserve(colvarcontent.size() + rcolvarcontent.size()-reducesize-1);//no append needed later
        mergedlist.insert(mergedlist.begin(), colvarcontent.begin(), colvarcontent.end());
        if (firstTime!=LastTime){ //if need to change the RC file
            for (std::size_t i=reducesize; i<rcolvarcontent.size(); ++i){
                std::string line = rcolvarcontent[i];
                std::stringstream stream2(line);
                double time;
                stream2>>time;
                std::string rest;
                std::getline(stream2,rest);
                std::ostringstream out;
                double corrected = time - firstTime + LastTime+PLUMEDSTEPSIZE;
                if (reducesize !=0){
                    out<<" ";
                }
                out<< corrected<<rest;
                mergedlist.push_back(out.str());
            }
        }else{
            mergedlist.insert(mergedlist.end(), rcolvarcontent.begin() +reducesize, rcolvarcontent.end());
        }
        std::filesystem::remove(newfolder +"/" + RCfile);
        std::filesystem::remove(newfolder+"/"+Cfile);
        std::ofstream combined(newfolder +"/"+ Cfile);
        for (const auto& line:mergedlist){
            combined<<line<<std::endl;
        }
    } 
}
int PLMD(std::string newfolder, std::string newfile, double leftbound, double rightbound, int runtype, int direction, int first, int strides){
    if (runtype<3){ //this is kind of neat, I'll put more detail in the manual on how this works, if basin does one handed run depending on case, if FFS does two handed run
        int elcaso = 3-2*runtype; //map 2 to -1 1 to 1
        if (elcaso*direction ==1){
            if (direction==1){
                leftbound = leftbound-100*(rightbound-leftbound); //unreachable
            }else{
                rightbound = rightbound + 100*(rightbound-leftbound);
            }
        }else{
            if (direction==1){
                rightbound = rightbound + 100*(rightbound-leftbound);
            }else{
                leftbound =leftbound-100*(rightbound-leftbound);
            }
        }
    }
    std::filesystem::copy_file("plumed.dat", newfolder + "/" + newfile +".dat", std::filesystem::copy_options::overwrite_existing);
    std::ifstream plumed(newfolder +"/" + newfile+".dat");
    std::vector<std::string> lines;//this was so much easier on python, the speedup better go crazy
    std::string line;
    int idxLL1 ;
    int idxUL1;
    int idxLL2;
    int idxUL2;
    int lineIndex=0 ;
    int hascommitor=0;
    while (std::getline(plumed,line)){
        if (line.find("FILE=") != std::string::npos) {
            std::stringstream stream(line);
            std::string token;
            std::string out;
            while (stream>>token){
                if (token.rfind("FILE=",0)==0){
                    if (first == 0) token = "FILE=" + std::filesystem::absolute(newfolder).string() + "/C" + newfile;
                    else token = "FILE=" + std::filesystem::absolute(newfolder).string() + "/RC" + newfile;
                }
                out += token + " ";
            }
            line=out;
        }
        if (line.find("STRIDE=") != std::string::npos){
            size_t indent = line.find_first_not_of(" \t");
            std::string prefix = (indent != std::string::npos) ? line.substr(0, indent) : "";

            std::stringstream stream(line);
            std::string token;
            std::string out;
            while (stream>>token){
                if (token.rfind("STRIDE=",0)==0){
                    token = "STRIDE=" + std::to_string(strides);
                }
                out += token + " ";
            }
            line = prefix + out;
        }
        if (line.find("BASIN_LL1=") != std::string::npos){
            idxLL1 =lineIndex;
        }
        if (line.find("BASIN_UL1=") != std::string::npos){
            idxUL1 =lineIndex;
        }
        if (line.find("BASIN_LL2=") != std::string::npos){
            idxLL2 =lineIndex;
        }
        if (line.find("BASIN_UL2=") != std::string::npos){
            idxUL2 =lineIndex;
        }
        lines.push_back(line);
        lineIndex++;
        if (line.find("COMMITTOR ...") != std::string::npos) { //write this after
            hascommitor+=1;
            if (runtype==4){
                lines.push_back("\tNOSTOP"); //don't stop on probe
                lineIndex++;
            }
        }
    }

    if (direction==1){
        lines[idxUL1]= changenumber(lines[idxUL1], leftbound);
        lines[idxLL2] = changenumber(lines[idxLL2], rightbound);
        lines[idxUL2] = changenumber(lines[idxUL2], rightbound +10*(rightbound-leftbound));
        lines[idxLL1] = changenumber(lines[idxLL1], leftbound-10*(rightbound-leftbound));
    }else{
        lines[idxUL1]= changenumber(lines[idxUL1], rightbound);
        lines[idxLL2] = changenumber(lines[idxLL2], leftbound);
        lines[idxUL2] = changenumber(lines[idxUL2], leftbound +10*(leftbound-rightbound));
        lines[idxLL1] = changenumber(lines[idxLL1], rightbound-10*(leftbound-rightbound));
    }
    plumed.close();
    std::ofstream PLUMED(newfolder +"/" + newfile+".dat");
    for (const auto&line:lines){
        PLUMED<<line<<std::endl;
    }
    PLUMED.close();
    //need to add runtypes
    //first copy the plumed.dat file into the folder with propername .dat
    return hascommitor;
}
double colvarruntime(std::string Colvar){
    std::string lastline = lastdataline(Colvar);
    std::stringstream stream(lastline);
    double time;
    stream>>time;
    return time;
}
std::pair<double, std::string> Runtimechecker(std::string newfolder){
    double runtime;
    std::string Path =newfolder+ "/path.txt";
    std::ifstream pathFile(Path);
    std::string propername ="";
    if (pathFile.is_open()){ 
        std::string Line;
        std::getline(pathFile, Line);
        std::vector<std::string> pathentries = pathlinereader(Line);
        propername = pathentries[2];
        if (std::filesystem::exists(newfolder + "/C" + propername)){
            runtime = colvarruntime(newfolder+"/C" + propername);
        }
        else{
            runtime =0;
        }

    }else{
        runtime=0;
    }
    return std::make_pair(runtime, propername); //apparently that's how you return two things in a function
}
void MDP(std::string newfolder, std::string propername, int seed, int calculator, double deltaT, long long nsteps ){
    std::filesystem::copy("ffs.mdp", newfolder +"/" +propername + ".mdp" , std::filesystem::copy_options::overwrite_existing);
    int velocitylinenumber = -1;

    //step 1 find the line that has velocity command
    std::string linetag;
    std::string deltatag;
    if (calculator ==0){
        linetag = "velocity";
        deltatag = "timestep";
    }else{
        linetag ="gen-seed";
        deltatag = "dt ";
    }
    std::ifstream mdp(newfolder +"/" + propername+".mdp");
    std::vector<std::string> lines;
    std::string line;
    int nstepsline = 0;
    int dtlinenumber=-1;
    int lineIndex=0 ; //borrowed this section from PLMD
    while (std::getline(mdp,line)){
        if (line.find(linetag) != std::string::npos){
            velocitylinenumber =lineIndex;
        }
        if (line.find(deltatag) !=std::string::npos){
            dtlinenumber = lineIndex;
        }
        if (line.find("nsteps ")!= std::string::npos){
            nstepsline = lineIndex;
        }
        lines.push_back(line);
        lineIndex++;
    } 
    mdp.close();
    if ((calculator ==0 && dtlinenumber!=-1 && velocitylinenumber==-1 ) || (calculator ==1 && dtlinenumber !=-1 && velocitylinenumber!=-1)) {
        if (calculator ==0){
            lines[dtlinenumber] = "timestep " + std::to_string(deltaT);
            std::ofstream mdp(newfolder +"/" + propername+".mdp");
            for (const auto&line:lines){
                mdp<<line<<std::endl;
            }
            mdp.close();
        }else{ //gromacs case
            size_t pos = lines[velocitylinenumber].find("= ");
            size_t posdt = lines[dtlinenumber].find("= ");
            size_t posnsteps = lines[nstepsline].find("= ");
            if (pos != std::string::npos && posdt != std::string::npos && posnsteps !=std::string::npos){
                lines[velocitylinenumber] = lines[velocitylinenumber].substr(0,pos) + "= " + std::to_string(seed);
                lines[dtlinenumber] = lines[dtlinenumber].substr(0,posdt) + "= " + std::to_string(deltaT) ;
                lines[nstepsline] = lines[nstepsline].substr(0,posnsteps) + "= " + std::to_string(nsteps);
                std::ofstream mdp(newfolder +"/" + propername+".mdp");
                for (const auto&line:lines){
                    mdp<<line<<std::endl;
                }
                mdp.close();
            }else{
                if (pos ==std::string::npos){
                    std::cout<<"Fix gen-seed line (add = sign)"<<std::endl;
                }if (posdt ==std::string::npos){
                    std::cout<<"Fix dt line (add = sign)"<<std::endl;
                }if (posnsteps ==std::string::npos){
                    std::cout<<"Fix nsteps line (add = sign)" <<std::endl;
                }
            }
        }
    }else{
        std::cout<<"Faulty MDP File"<<std::endl;
    }
}
void MDPrestart(std::string newfolder, std::string propername, int calculator){
    //normally I'd like to handle the similar parts together but in this one there's a different number of operations so easier to split
    if (calculator ==0){
        int velocitylinenumber = -1;
        std::ifstream mdp(newfolder +"/" + propername+".mdp"); //had this from when I was using velocity initialization, but it shouldn't fire either way so it's okay to leave it
        std::vector<std::string> lines;
        std::string line;
        int lineIndex=0 ; //borrowed this section from PLMD
        while (std::getline(mdp,line)){
            if (line.find("velocity") != std::string::npos){
                velocitylinenumber =lineIndex;
            }
            lines.push_back(line);
            lineIndex++;
        }
        if (velocitylinenumber!=-1){//if there is no line don't need to do anything
            line = lines[velocitylinenumber];
            int pos = line.find_first_not_of(" \t");
            if (line[pos] != '#'){ //if commented don't need to do anything
                lines[velocitylinenumber] = "#" + line ;
                mdp.close();
                std::ofstream mdp(newfolder +"/" + propername+".mdp");
                for (const auto&line:lines){
                    mdp<<line<<std::endl;
                }
                mdp.close();
            }
        }
    }else{
        int line1index=-1;
        int line2index = -1;
        int line3index =-1;
        std::ifstream mdp(newfolder +"/" + propername + ".mdp");
        std::vector<std::string> lines;
        std::string line;
        int lineIndex = 0;
        while (std::getline(mdp, line)){
            if (line.find("gen-vel") != std::string::npos){
                line1index = lineIndex;
            }
            if (line.find("gen-temp")!= std::string::npos){
                line2index = lineIndex;
            }
            if (line.find("gen-seed")!= std::string::npos){
                line3index = lineIndex;
            }
            lines.push_back(line);
            lineIndex++;
        }
        size_t pos2 = lines[line2index].find_first_not_of(" \t");
        if (pos2 != std::string::npos && lines[line2index][pos2] != ';'){
            lines[line2index] = ";" + lines[line2index];
        }
        size_t pos3 = lines[line3index].find_first_not_of(" \t");
        if (pos3 != std::string::npos && lines[line3index][pos3] != ';'){
            lines[line3index] = ";" + lines[line3index];
        }
        size_t pos = lines[line1index].find("= ");
        if (pos != std::string::npos){
            lines[line1index] = lines[line1index].substr(0,pos) + "= no";
            std::ofstream mdp(newfolder +"/" + propername+".mdp");
            for (const auto&line:lines){
                mdp<<line<<std::endl;
            }
            mdp.close();
        }else{
            std::cout<<"Fix gen-seed line (add = sign)"<<std::endl;
        }   
    }
}
void deleteall(const std::string &folder)
{
    for (const auto &entry : std::filesystem::directory_iterator(folder)) {
        std::filesystem::remove_all(entry.path());
    }
}

void wipe(std::string folder, std::string propername) //remove those files in the subdirectory that are from a certain path line
{
    for (const auto &entry : std::filesystem::directory_iterator(folder)) {
        std::string name = entry.path().filename().string();
        if (name.find(propername) != std::string::npos) {
            std::filesystem::remove_all(entry.path());
        }
    }
}



int Repeater(std::string newfolder, int index, int total){
    std::ifstream inputFile(newfolder +"/iterationtracker.txt");
    int number=0;
    if (!inputFile.is_open()){
        //if file doesn't exist populate it entirely with zeros and one 1 at index
        std::ofstream OutFile(newfolder + "/iterationtracker.txt");
        for (int i=0; i<total; i++){
            if (i ==index){
                OutFile<< 1<<std::endl;
            }
            else{
                OutFile<< 0<<std::endl;
            }
        }
        number =1;
        OutFile.close(); 
    }else{
        int Prior;
        std::vector<std::string> lines;
        for (int i=0; i<total; i++){
            inputFile>>Prior;
            if (i ==index){
                number = Prior +1;
                lines.push_back( std::to_string(Prior +1) +"\n");
            }else{
                lines.push_back(std::to_string(Prior) +"\n");
            }
        }
        inputFile.close();
        std::ofstream OutFile(newfolder + "/iterationtracker.txt");
        for (int i =0; i<total; i++){
            OutFile<<lines[i];
        }
    }
    return number;    
}
int gronameparser(std::string structurefiles, int elcual){
    std::ifstream groname(structurefiles + "/gronames.txt");
    std::string Line;
    int out;
    int tag;
    std::string ic = "initialconfiguration_";
    while (std::getline(groname, Line)){
        size_t pos = Line.find_last_of(',');
        if (std::stoi(Line.substr(pos +1))==elcual){
            size_t start = Line.find(ic);
            start += ic.length();
            size_t end = Line.find("_",start);
            out = std::stoi(Line.substr(start, end-start));
        }
    }
    return out;
}
LAMMPS_NS::LAMMPS* firstgo(Context* ctx , std::string structurefiles, std::string Outputfolder, std::string runid, double leftbound, double rightbound, int runtype, int direction,  std::string CUSTOMOP, OrderParameter* op, int iteration, long long nsteps , std::string trail,std::string rootdirectory, int calculator, int strides, double deltaT, int constanttopology){
    const int rank = ctx->rank();
    int hascommitor;
    std::hash<std::string> hasher;//different for each run
    int seed = static_cast<int>(hasher("folder" + Outputfolder + "slurm" + runid + "iteration"+std::to_string(iteration)) & 0x7FFFFFFF); // in FFS must specify iteration
    std::mt19937 rng(seed); //alright now we have reproducibility
    int usePLMD = 0; 
    if (CUSTOMOP.empty()){
        usePLMD = 1;
    }
    std::string newfolder = Outputfolder +"/" + runid ;
    std::vector<std::string> inputfiles;
    std::string inputfile;//in case somehow the pick is different
    int repeat = 1; //must have FFS tell which repeat this is
    if (rank ==0){
        inputfiles = listdirEnding(structurefiles, trail);
        if (runtype==3){//probabilistic
            std::uniform_int_distribution<size_t> dist(0, inputfiles.size() - 1);
            int whichindex = dist(rng);
            repeat = Repeater(newfolder, whichindex, inputfiles.size());
            inputfile = inputfiles[whichindex];
        }else{//basin or exploration
            inputfile = inputfiles[std::stoi(runid)-1];
        }
    }
    ctx->barrier();
    ctx->bcast_string(inputfile, 0);
    ctx->bcast_int(repeat, 0);


    size_t pos = inputfile.rfind('.');
    std::string beforedot;
    beforedot = inputfile.substr(0,pos);
    std::string propername = "folder_" + Outputfolder + "_subfolder_" + runid +"_initialconfiguration_" + beforedot +"_iteration_"+std::to_string(repeat) ;//must be one in basin
    
    //write all inputfiles into a txt file called grofiles
    if (rank ==0){
        if (!std::filesystem::exists(newfolder + "/possibleinputs.txt")){
            std::ofstream Ins(newfolder +"/possibleinputs.txt");
            for (const std::string& entry:inputfiles){
                Ins<< entry<< std::endl;
            }
            Ins.close();
        }
        //write chosen one to outputfiles
        std::ofstream Used(newfolder +"/outputfiles.txt", std::ios::app); //this also has to be an appending
        Used<< inputfile<<std::endl;
        Used.close();

        //no need for iteration tracker in basin code
        //set random seed and save it
        std::ofstream SeedTXT(newfolder +"/seedvalues.txt", std::ios::app); //thjis as well
        SeedTXT <<std::to_string(seed) <<std::endl;
        SeedTXT.close();
        //make a path file
        std::ofstream PathFile(newfolder +"/path.txt", std::ios::app);//have to append to it with each one
        PathFile<< beforedot <<", " << std::to_string(repeat)<<", "<<propername<<std::endl; //in python I woulda done this with like str + str I think in c++ it's <<
        //assume run parameters are kept in an ffs.mdp style file
        MDP(newfolder, propername, seed, calculator,deltaT, nsteps); //this should have all kinds of mdp file edits but will have to see what these look like in LAAMPS implementation first.
        //take plumed file modify output locations and modify the order parameter values
        if (usePLMD ==1){
            hascommitor=PLMD(newfolder, propername, leftbound, rightbound,runtype,direction,0, strides);
        }
    }
    ctx->barrier();
    ctx->bcast_int(hascommitor, 0);
    if (usePLMD==1){
        if (hascommitor==0){
            std::cout<<"plumed missing commitor"<<std::endl;
            return 0;
        }
    }
    LAMMPS_NS::LAMMPS* lmp = nullptr;
    //randomize seed, come back to this 
    //Here for gromacs I delete all the unncessary files come back to here if there are any in LAMMPS
    if (calculator ==0){
#ifdef USE_LAMMPS
        lmp = InitializeLammps(ctx,newfolder, propername, structurefiles+"/"+inputfile, usePLMD); //gotta add the directory here as I removed it from listdir
        //gotta get number of steps as an input, later will check colvar files to see how many steps are already completed.
        lammpsrunner(ctx,lmp, op, usePLMD, nsteps, 0, direction, runtype,leftbound, rightbound, newfolder, propername, strides,deltaT );
        writeCPT(lmp, newfolder,propername);
        ctx->barrier();
#endif
    }else{
        if (trail.find("gro") != std::string::npos){
            //if this is the case we don't have to feed in a cpt file
            std::string root = newfolder + "/" + propername;
            std::string substring = inputfile.substr(0, inputfile.size() - 4);
            if (constanttopology ==0){
                InitializeGMX(root +  ".mdp",structurefiles +"/" + inputfile, structurefiles + "/" + substring + ".top",structurefiles + "/"+ substring+ ".ndx", root, "",ctx);
            }else{
                InitializeGMX(root +  ".mdp",structurefiles +"/" + inputfile, "topol.top","index.ndx", root, "",ctx);
            }
            ctx->barrier();
        }else{ //I'm going to have the safety aspects be in the preliminary user inputs scripts, this way I can assume there's a gronames.txt file here etc
            //step 1 open gronames.txt 2. find line with correct input file, 3. split to find first _initialconfiguration_X_, 4. Feed X.gro with path as .gro, 5. feed input as cpt 
            std::string root = newfolder + "/" + propername;
            int grofile = gronameparser(structurefiles, std::stoi(beforedot));
            if (constanttopology==0){
                InitializeGMX(root + ".mdp",rootdirectory + "/"+std::to_string(grofile) +".gro",  rootdirectory + "/"+std::to_string(grofile)+".top", rootdirectory + "/"+std::to_string(grofile) + ".ndx", root, structurefiles + "/" + inputfile,ctx);
            }else{
                InitializeGMX(root + ".mdp",rootdirectory + "/"+std::to_string(grofile) +".gro",  "topol.top", "index.ndx", root, structurefiles + "/" + inputfile,ctx);
            }
            ctx->barrier();
        }
        //left off here, now run the gmx sim
        RUNGMX(newfolder + "/" + propername,0,ctx);
        ctx->barrier();
    }
    return lmp;
}


void RecordCrossing(std::string folder, double time, int crossings){
    int linecount= 0;
    std::cout<<"Number of Crossings: "<<crossings<<std::endl;
    if (std::filesystem::exists(folder + "/BasinCrossingRecorder.txt")){
        std::ifstream crossingrecord(folder+ "/BasinCrossingRecorder.txt");
        std::string line;
        while (std::getline(crossingrecord, line)){
            linecount++; //I wanted to make sure the user has some file that gives when each crossing happened
        }// in the previous python version this was build in to basininfo.txt but it felt kinda odd having them both together
        crossingrecord.close();
    }
    
    if (linecount<crossings){
        std::ofstream crossingrecord(folder + "/BasinCrossingRecorder.txt", std::ios::app);
        crossingrecord<<"Crossing: "<< crossings <<", "<<"Time: "<<time<<std::endl;
    }

}

void pathedit(std::string newfolder,int& runtype, int direction, double leftbound, double rightbound,int& weshouldstop, long long nsteps, double deltaT, int targetline ){
    std::vector<std::string> pathLines;
    std::string propername;
    std::vector<std::string> pathentries;
    std::string Path = newfolder + "/path.txt";
    std::ifstream pathFile(Path);
    double c;
    std::string Line;
    int previouscrossings;
    int currentline=-1;
    while (std::getline(pathFile, Line)){
        currentline++;
        if (runtype==3){ //ffs operation
            pathentries = pathlinereader(Line);
            if (pathentries.back().find("success")!= std::string::npos || pathentries.back().find("failure")!= std::string::npos){
                pathLines.push_back(Line); //was already marked as a finished simulation, no need to check further
                if (currentline==targetline){
                    weshouldstop = 0;
                }
            }else{
                propername =pathentries[2];
                std::string Colvar = newfolder + "/C" + propername;
                if (std::filesystem::exists(Colvar)){
                    auto result = lastorderparameter(Colvar);
                    double time = result.first; //no need to convert to nsteps here
                    c= result.second;
                    if (c*direction>rightbound*direction){
                        pathLines.push_back( pathentries[0] + ", "+pathentries[1] + ", " + pathentries[2] + ", success");
                        if (currentline==targetline){
                            weshouldstop = 0;
                        }
                    }else if (c*direction<leftbound*direction){
                        pathLines.push_back( pathentries[0] + ", "+pathentries[1] + ", " + pathentries[2] + ", failure");
                        if (currentline==targetline){
                            weshouldstop = 0;
                        }
                    }else{
                        pathLines.push_back( pathentries[0] + ", "+pathentries[1] + ", " + pathentries[2] + ", incomplete");
                    }
                }else{
                    wipe(newfolder, propername);  //problem case, delete colvar and delete from path
                    std::cout<<"colvar not developing"<<std::endl;
                    weshouldstop =2; //run issue
                }
            }
        }else if(runtype==2 || runtype ==1){//there's actually no need to check if the colvar exists as we wouldn't have entered otherwise
            pathentries = pathlinereader(Line);
            propername =pathentries[2];
            std::string Colvar = newfolder + "/C" + propername;
            if (pathentries.size()<4){
                previouscrossings=0;
            }else{
                previouscrossings = std::stoi(pathentries[3]);
            }
            int nocrossings = crossings_checker(Colvar, leftbound,rightbound, direction); 
            if (nocrossings>previouscrossings){
                weshouldstop =0;
            }
            pathLines.push_back( pathentries[0] + ", " + pathentries[1] + ", " + pathentries[2] + ", " + std::to_string(nocrossings));
            
            std::string basindirectoryfile = newfolder +"/basininfo.txt";
            auto result = lastorderparameter(Colvar);
            double time = result.first; //no need to convert to nsteps here
            c= result.second;
            if (direction*c >= direction *rightbound && runtype ==1){
                std::filesystem::copy(newfolder + "/" + propername +".cpt", newfolder + "/" + std::to_string(nocrossings) + ".cpt",std::filesystem::copy_options::overwrite_existing); //export crossing
                RecordCrossing(newfolder, time,nocrossings);
                runtype =2;
                weshouldstop = 0;
            } else if (direction*c <= direction*leftbound && runtype ==2){
                weshouldstop = 0;
                runtype = 1;
                //switch case
                //imo only have to switch case and edit basininfo but why not also edit the path to the right number of crossings just in case
            }
        std::ofstream basinOut(basindirectoryfile, std::ios::trunc);
        basinOut<<runtype<<std::endl;
        }else if (runtype ==4){ //this is for the exploration
            pathentries = pathlinereader(Line);
            auto result = Runtimechecker(newfolder);
            double runtime   = result.first/deltaT;
            std::string propername = result.second;
            if (runtime>=nsteps){//will need to wire in nsteps
                pathLines.push_back( pathentries[0] + ", " + pathentries[1] + ", " + pathentries[2] + ", complete");
                weshouldstop = 0;
            }
            else{
                pathLines.push_back(Line);
            }
        }
    }
    std::ofstream pathOut(Path, std::ios::trunc);
    for (std::string line:pathLines){
        pathOut<<line<<std::endl;
    }
}

void DeleteExtraFiles(std::string newfolder, int runtype){
    if (runtype ==3){//ffs operation
        std::vector<std::string> HashFiles = listdirStart(newfolder, "#");
        for (const auto&file : HashFiles){
            std::filesystem::path p = std::filesystem::path(newfolder) /file;
            if (std::filesystem::exists(p)){
                std::filesystem::remove(p);
            }
        }
        std::vector<std::string> endings = {".edr", ".log", ".tpr", ".trr"};
        for (const auto&ending: endings){
            std::vector<std::string> TRRs = listdirEnding(newfolder, ending);
            for (const auto& file: TRRs){
                std::filesystem::path p = std::filesystem::path(newfolder) / file;
                if (std::filesystem::exists(p)){
                    std::filesystem::remove(p);
                }
            }
        }
    }
    std::vector<std::string> starts = {"bck."}; //leaving it open ended in case have to add more
    for (const auto& start:starts){
        std::vector<std::string> BCKS = listdirStart(newfolder, start);
        for (const auto& file: BCKS){
            std::filesystem::path p = std::filesystem::path(newfolder) / file;
            if (std::filesystem::exists(p)){
                std::filesystem::remove(p);
                std::cout<<"Remove "<<p<<std::endl;
            }
        } 
    }
}

void CleanUp(std::string newfolder, int& runtype, int direction, double leftbound, double rightbound,int& weshouldstop, long long nsteps, double strides, double deltaT, int calculator , int usePLMD, int targetline){
    double PLMDSTEPSIZE = strides*deltaT;
    DeleteExtraFiles(newfolder, runtype);
    int feedintoRC = 2*usePLMD-calculator; //corelary to how many lines should skip when appending RC, gromacs with PLMD skip 1 line, lammps with plumed skip 2, lammps without skip nothing
    std::vector<std::string> restartcolvars = listdirStart(newfolder, "RC");//RC file combiner, could do this via last entry and path but this is clean enough and will translate better in FFS
    if (!restartcolvars.empty()){
        for (std::string RC : restartcolvars){
            RCfilecombiner(newfolder, RC,PLMDSTEPSIZE,feedintoRC); //append RC file to existing C file
        }
    }
    std::vector<std::string> colvars = listdirStart(newfolder, "C");
    weshouldstop = 1;
    if (colvars.size()==0){
        deleteall(newfolder); //there's some issue if no colvar
        weshouldstop =2; //this will cause termination making it so the script doesn't requeue
    }else{
        pathedit(newfolder, runtype, direction, leftbound, rightbound, weshouldstop, nsteps,deltaT, targetline);
    }
}

int ProbabilisticCompletion(std::string newfolder){
    int count = 0;
    std::string pathFile = newfolder + "/path.txt";
    if (std::filesystem::exists(pathFile)){
        std::ifstream Path(pathFile);
        std::string Line;
        while (std::getline(Path, Line)){
            if (Line.find("success")!= std::string::npos || Line.find("failure")!= std::string::npos){
                count ++;
            }
        }
        Path.close();
    }else{
        return -1;
    }
    return count;
}
void TerminationStall(int var){
    if (var ==2){
        std::cout<< "Run Initiation Failure"<<std::endl; //note to self create error codes
        std::exit(EXIT_FAILURE); //in case colvars aren't developing
    }
    if (var==1){
        std::cout<<"Sleeping"<<std::endl;
    }
    while (var == 1) {
        std::this_thread::sleep_for(std::chrono::seconds(60));
    }
}

void Basin(Context* ctx,std::string Outputfolder, std::string runid, double leftbound, double rightbound, std::string structurefiles, std::string trail, long long nsteps, std::string CUSTOMOP, std::string rootdirectory, int calculator, int strides, double deltaT, int constanttopology){
    LAMMPS_NS::LAMMPS* lmp = nullptr;
    int hascommitor;
    int usePLMD = 0;
    if (CUSTOMOP.empty()){
        usePLMD = 1;
    }
    const int rank = ctx->rank();
    std::unique_ptr<OrderParameter> op; //check later on if it's okay for this to live outside the scope
    if (usePLMD ==0){
        op = createOrderParameter(CUSTOMOP); //alright this should let me take an input string and use it to map to an op calculating function
    }
    //make non principle threads idle until lammps step


    // let's break down the proccess of running a basin code into steps, I have a pretty good workflow in gromacs will mimc
    // the first step naturally is to create output folders and so on
    //the second step will be to define the direction of the basin run like a half FFS step etc.
    // then we enter the action loop
    //step 1 make folders, get direction of basin run
    std::string newfolder = Outputfolder +"/" + runid ;
    std::string basindirectoryfile = newfolder + "/basininfo.txt";
    int runtype=1;
    if (!std::filesystem::exists(newfolder)){
        if (rank ==0){
            std::filesystem::create_directories(newfolder);
        }
        ctx->barrier();
        runtype = 1;
    }else{
        if (rank ==0){
            std::ifstream inputFile(basindirectoryfile);
            if (!inputFile.is_open()){
                runtype =1; //if file can't open for any reason default to standard runtype, including if file doesn't exist
            }else{
                inputFile>> runtype;
                inputFile.close();
            }
        }
        ctx->barrier();
        ctx->bcast_int(runtype, 0);
    }
    int direction =rightbound>leftbound?1:-1; //pos 1 for increasing, neg one for decreasing
    
    std::string Path =newfolder+ "/path.txt"; //preemptive edits to files in case of interruptions
    int trashvariable = 0; //idk colaboratory with cleanup
    if (std::filesystem::exists(Path)){
        if (rank ==0){
            CleanUp(newfolder, runtype, direction, leftbound, rightbound, trashvariable,nsteps, strides,deltaT,calculator, usePLMD,0);
        }
    }
    ctx->barrier(); 
    ctx->bcast_int(runtype, 0);
    //check if is first launch
    int pathlines = 0;
    std::string Line;
    if (rank ==0){
        std::ifstream pathfile(Path);
        if (pathfile.is_open()){
            while (std::getline(pathfile,Line)){
                pathlines++;
            }
            pathfile.close();
        }
    }
    ctx->barrier();
    ctx->bcast_int(pathlines, 0);    
    int weshouldstop = 0;
    if (pathlines==0){//okay now the script for first launch basin
        lmp=firstgo(ctx,structurefiles, Outputfolder, runid, leftbound, rightbound, runtype, direction,CUSTOMOP,op.get(),1,nsteps,trail, rootdirectory, calculator, strides,deltaT,constanttopology); //so I know it's kind of weird to have this be a function but it's because I need to do this twice, saves on number of lines
        //now clean the files
        if (rank ==0){
            CleanUp(newfolder, runtype, direction, leftbound, rightbound, weshouldstop,nsteps, strides,deltaT,calculator, usePLMD,0);
        }
        ctx->barrier();
        ctx->bcast_int(weshouldstop, 0);
        ctx->bcast_int(runtype, 0);
        DestroyLammps(lmp);
        TerminationStall(weshouldstop); //so when a termination happens the script keeps running a little while but gmx and so on commands get halted, to avoid issues it's best to check if term happened and stop making attempts if so.
    }
    long long runtime;
    std::string propername;
    if (rank ==0){
        auto result = Runtimechecker(newfolder);
        runtime   = result.first/deltaT;
        propername = result.second;
    }
    ctx->barrier();
    ctx->bcast_longlong(runtime, 0);
    ctx->bcast_string(propername, 0);
    weshouldstop = 0;
    while (runtime <nsteps){
        //asssume already have correct runtime,propername  (calculate before while and at end of loop)
        //step one, trigger a restart run
        //this involves removing mdp vel randomization, changing plumed, deleting unneeded files, and triggering a run command
        //step 1 mdp randomization on LAMMPS
        if (std::filesystem::exists(newfolder +"/C" + propername) and std::filesystem::exists(newfolder + "/" + propername + ".cpt")){
            if (rank ==0){
                MDPrestart(newfolder, propername,calculator);
                //step 2 plumed editing 
                if (usePLMD ==1){
                    hascommitor= PLMD(newfolder, propername,leftbound,rightbound,runtype, direction,1, strides);
                }
            }
            ctx->barrier();
            ctx->bcast_int(hascommitor, 0);
            if (usePLMD ==1){
                if (hascommitor ==0){
                    std::cout<<"plumed missing commitor"<<std::endl;
                    return;
                }
            }
            if (calculator ==0){
#ifdef USE_LAMMPS
                lmp = RestartLammps(ctx,newfolder, propername, usePLMD);
                lammpsrunner(ctx, lmp, op.get(), usePLMD, nsteps, runtime, direction, runtype, leftbound, rightbound, newfolder, propername, strides,deltaT);
                writeCPT(lmp, newfolder, propername);
                ctx->barrier();
#endif
            }else{
                RUNGMX(newfolder + "/" + propername, 1,ctx);
                ctx->barrier();
            }
        
        }else{
            //if checkpoint doesn't exist we should delete colvars and then start from scratch, also reset steps etc
            //step 1, reset path file, basininfo etc, step 2, start a new. I think what I'll do is write firsttry as a function and have this and the previous do it
            //will just add here an extra delete all files in folder command
            if (rank ==0){
                deleteall(newfolder);
            }
            ctx->barrier();
            runtime = 0;
            runtype =1;
            lmp=firstgo(ctx,structurefiles, Outputfolder, runid, leftbound, rightbound, runtype, direction,CUSTOMOP,op.get(),1,nsteps,trail, rootdirectory, calculator, strides,deltaT,constanttopology);
        }
        if (rank ==0){
            CleanUp(newfolder, runtype, direction, leftbound, rightbound, weshouldstop,nsteps, strides,deltaT,calculator, usePLMD,0);
            auto result = Runtimechecker(newfolder);
            runtime   = result.first/deltaT;
            propername = result.second;
        }//this is a rather long one make sure all things here are fine, must have runtime and propername transfer and have return statement if colvar doesn't develop
        ctx->barrier();
        ctx->bcast_longlong(runtime, 0);
        ctx->bcast_string(propername, 0);
        ctx->bcast_int(weshouldstop, 0);
        ctx->bcast_int(runtype, 0);
        DestroyLammps(lmp);
        if (runtime<nsteps){
            TerminationStall(weshouldstop); //so another way we can terminate without plumed intervention causing weshouldstop to trigger termination stall is if we hit nsteps, this is a failsafe to allow the code to exit
        }
    }
}


void Probabilistic(Context* ctx,std::string Outputfolder, std::string runid, double leftbound, double rightbound, int desiredtrajectories, std::string structurefiles, std::string trail, std::string CUSTOMOP, std::string rootdirectory, int calculator, int strides, double deltaT, int constanttopology){
    LAMMPS_NS::LAMMPS* lmp= nullptr;
    std::vector<std::string> pathentries;
    int hascommitor;
    long long nsteps =std::numeric_limits<int>::max();
    int usePLMD = 0;
    int runtype =3;
    if (CUSTOMOP.empty()){
        usePLMD = 1;
    }
    const int rank = ctx->rank();
    std::unique_ptr<OrderParameter> op;
    if (usePLMD ==0){
        op = createOrderParameter(CUSTOMOP) ;//alright this should let me take an input string and use it to map to an op calculating function
    }
    //make non principle threads idle until lammps step

    std::string newfolder = Outputfolder +"/" + runid ;
    if (!std::filesystem::exists(newfolder)){
        if (rank ==0){
            std::filesystem::create_directories(newfolder);
        }
        ctx->barrier();
    }
    int direction =rightbound>leftbound?1:-1; //pos 1 for increasing, neg one for decreasing
    std::string Path =newfolder+ "/path.txt"; //preemptive edits to files in case of interruptions
    int trashvariable = 0;
    if (std::filesystem::exists(Path)){
        if (rank ==0){
            CleanUp(newfolder, runtype, direction, leftbound, rightbound, trashvariable,nsteps, strides,deltaT,calculator, usePLMD,0);
        }
    }
    ctx->barrier();
    
    //check if is first launch for any
    int pathlines = 0;
    std::string Line;
    if (rank ==0){
        std::ifstream pathfile(Path);
        if (pathfile.is_open()){
            while (std::getline(pathfile,Line)){
                pathlines++;
            }
            pathfile.close();
        }
    }
    ctx->barrier();
    ctx->bcast_int(pathlines, 0);
    int weshouldstop = 0; 
    for (int iteration=pathlines; iteration<desiredtrajectories; iteration++){
        lmp=firstgo(ctx,structurefiles, Outputfolder, runid, leftbound, rightbound, 3 , direction, CUSTOMOP, op.get(),iteration,nsteps,trail, rootdirectory, calculator, strides,deltaT,constanttopology); //this will attempt to complete FFS from first launch for each of the desired trajectories
        if (rank ==0){
            CleanUp(newfolder, runtype, direction, leftbound, rightbound, weshouldstop,nsteps, strides,deltaT,calculator, usePLMD, iteration);
        }
        ctx->barrier();
        ctx->bcast_int(weshouldstop, 0); //if it weren't for the scavenge partition we could probably leave it at this
        DestroyLammps(lmp);
        TerminationStall(weshouldstop);
    }

    int Completed=0;
    if (rank ==0){
        Completed = ProbabilisticCompletion(newfolder);
    }
    ctx->barrier();
    ctx->bcast_int(Completed, 0);

    if (Completed==-1){
        std::cout<<"Path Didn't Generate"<<std::endl;
        return;
    }
    std::string PathLine;
    int attempts=0;
    std::string propername;
    std::ifstream readpath(Path);
    if (Completed!=desiredtrajectories){
        for (int trajectory=0; trajectory<desiredtrajectories; trajectory++){
            if (std::getline(readpath, PathLine)){
                pathentries = pathlinereader(PathLine);
                if (pathentries.back().find("success")== std::string::npos && pathentries.back().find("failure")== std::string::npos){
                    propername = pathentries[2];
                    if (rank ==0){
                        MDPrestart(newfolder, propername,calculator);
                        if (usePLMD ==1){
                            hascommitor=PLMD(newfolder, propername,leftbound, rightbound, runtype, direction,1, strides);
                        }
                    }
                    ctx->barrier();
                    ctx->bcast_int(hascommitor, 0);
                    if (usePLMD==1){
                        if (hascommitor==0){
                            std::cout<<"plumed missing commitor"<<std::endl;
                            return;
                        }
                    }
                    if (calculator ==0){
#ifdef USE_LAMMPS
                        lmp = RestartLammps(ctx,newfolder, propername, usePLMD);
                        lammpsrunner(ctx,lmp,op.get(), usePLMD, nsteps, 0, direction, runtype, leftbound, rightbound, newfolder, propername, strides,deltaT);
                        writeCPT(lmp, newfolder, propername);
                        ctx->barrier();
#endif
                    }else{
                        std::string root = newfolder + "/" + propername;
                        std::string ic = "initialconfiguration_";
                        size_t start =propername.find(ic); //need to extract ic from propername
                        start += ic.length();
                        size_t end = propername.find("_", start);
                        int chosenstructure = std::stoi(propername.substr(start,end-start));

                        if (trail.find("gro")!=std::string::npos){
                            if (constanttopology ==0){
                                InitializeGMX(root + ".mdp", structurefiles +"/" + std::to_string(chosenstructure) + ".gro", structurefiles + "/" + std::to_string(chosenstructure) + ".top", structurefiles +"/" + std::to_string(chosenstructure) + ".ndx", root, root + ".cpt",ctx);   
                            }
                            else{
                                InitializeGMX(root + ".mdp", structurefiles +"/" + std::to_string(chosenstructure) + ".gro",  "topol.top",  "index.ndx", root, root + ".cpt",ctx);
                            }
                            ctx->barrier();
                        }else{ //I'm going to have the safety aspects be in the preliminary user inputs scripts, this way I can assume there's a gronames.txt file here etc
                            //step 1 open gronames.txt 2. find line with correct input file, 3. split to find first _initialconfiguration_X_, 4. Feed X.gro with path as .gro, 5. feed input as cpt 
                            int grofile = gronameparser(structurefiles, chosenstructure);
                            if (constanttopology==0){
                                InitializeGMX(root + ".mdp",rootdirectory + "/"+std::to_string(grofile) +".gro",  rootdirectory +"/"+ std::to_string(grofile)+".top", rootdirectory + "/"+std::to_string(grofile) + ".ndx", root, root + ".cpt",ctx);
                            }else{
                                InitializeGMX(root + ".mdp",rootdirectory + "/"+std::to_string(grofile) +".gro", "topol.top", "index.ndx", root, root + ".cpt",ctx);
                            }
                            ctx->barrier();
                        }
                        RUNGMX(root, 0, ctx);
                    }
                }
            }else{//with any luck this should never fire
                lmp=firstgo(ctx,structurefiles, Outputfolder, runid, leftbound, rightbound, 3 , direction, CUSTOMOP, op.get(),trajectory,nsteps,trail, rootdirectory, calculator, strides,deltaT,constanttopology);
            }
            //now clean
            if (rank ==0){
                CleanUp(newfolder, runtype, direction, leftbound, rightbound, weshouldstop,nsteps, strides,deltaT,calculator, usePLMD,trajectory);
            }
            ctx->barrier();
            ctx->bcast_int(weshouldstop, 0);
            DestroyLammps(lmp);
            TerminationStall(weshouldstop);
        }
    }
    //here i'm assuming if timeout makes exit loop it will not then count as complete
}


void Exploration(Context* ctx, std::string Outputfolder, std::string runid, std::string structurefiles, std::string trail, long long nsteps, std::string CUSTOMOP, std::string rootdirectory, int calculator, int strides, double deltaT, int constanttopology){
    LAMMPS_NS::LAMMPS* lmp=nullptr;
    int hascommitor;
    int runtype =4;
    int usePLMD = 0;
    double rightbound = 10;
    double leftbound = 1;
    int direction = 1; //fill these in w lb, rb placeholders as they wont matter
    if (CUSTOMOP.empty()){
        usePLMD = 1;
    }
    //parallel context ranks
    const int rank = ctx->rank();

    std::unique_ptr<OrderParameter> op;
    if (usePLMD ==0){
        op = createOrderParameter(CUSTOMOP); //alright this should let me take an input string and use it to map to an op calculating function
    }
    //make non principle threads idle until lammps step

    //this is going to be basically the basin exploration but with no lb, rb
    std::string newfolder = Outputfolder +"/" + runid ;
    if (!std::filesystem::exists(newfolder)){
        if (rank ==0){
            std::filesystem::create_directories(newfolder);
        }
        ctx->barrier();
    }    
    std::string Path =newfolder+ "/path.txt"; //preemptive edits to files in case of interruptions
    int trashvariable = 0; //idk colaboratory with cleanup
    if (std::filesystem::exists(Path)){
        if (rank ==0){
            CleanUp(newfolder, runtype, direction, leftbound, rightbound, trashvariable,nsteps, strides,deltaT,calculator, usePLMD,0); 
        }
    }
    ctx->barrier(); //I dont think anything must be communicated here
    //check if is first launch
    int pathlines = 0;
    std::string Line;
    if (rank ==0){
        std::ifstream pathfile(Path);
        if (pathfile.is_open()){
            while (std::getline(pathfile,Line)){
                pathlines++;
            }
            pathfile.close();
        }
    }
    ctx->barrier();
    ctx->bcast_int(pathlines, 0);
    int weshouldstop =0;
    if (pathlines==0){//okay now the script for first launch basin
        lmp=firstgo(ctx,structurefiles, Outputfolder, runid, leftbound, rightbound, runtype, direction,CUSTOMOP,op.get(),1,nsteps,trail, rootdirectory, calculator, strides,deltaT,constanttopology); //so I know it's kind of weird to have this be a function but it's because I need to do this twice, saves on number of lines
        //now clean the files
        if (rank ==0){
            CleanUp(newfolder, runtype, direction, leftbound, rightbound, weshouldstop,nsteps, strides,deltaT,calculator, usePLMD,0);
        }
        ctx->barrier();
        ctx->bcast_int(weshouldstop, 0);
        DestroyLammps(lmp);
        TerminationStall(weshouldstop);
    }
    long long runtime;
    std::string propername;
    if (rank ==0){
        auto result = Runtimechecker(newfolder);
        runtime   = result.first/deltaT;
        propername = result.second;
    }
    ctx->barrier();
    ctx->bcast_longlong(runtime, 0);
    ctx->bcast_string(propername, 0);

    while (runtime <nsteps){
        if (std::filesystem::exists(newfolder +"/C" + propername) && std::filesystem::exists(newfolder + "/"+propername +".cpt")){
            //step 1 mdp randomization on LAMMPS
            if (rank ==0){
                MDPrestart(newfolder, propername, calculator);
                //step 2 plumed editing 
                if (usePLMD ==1){
                    hascommitor=PLMD(newfolder, propername,leftbound,rightbound,runtype, direction,1, strides);
                }
                //step 3 run command
            }
            ctx->barrier();
            ctx->bcast_int(hascommitor, 0);
            if (usePLMD==1){
                if (hascommitor==0){
                    std::cout<<"plumed missing commitor"<<std::endl;
                    return;
                }
            }
            if (calculator ==0){
#ifdef USE_LAMMPS
                lmp = RestartLammps(ctx,newfolder, propername,usePLMD); // no gromacs grompp command in this case
                lammpsrunner(ctx,lmp, op.get(), usePLMD, nsteps, runtime, direction, runtype, leftbound, rightbound, newfolder, propername, strides,deltaT);//now run remaining steps
                writeCPT(lmp, newfolder, propername);
#endif
            }else{
                RUNGMX(newfolder + "/" + propername,1,ctx);
                ctx->barrier();
            }
        
        }else{
            //if checkpoint doesn't exist we should delete colvars and then start from scratch, also reset steps etc
            //step 1, reset path file, basininfo etc, step 2, start a new. I think what I'll do is write firsttry as a function and have this and the previous do it
            //will just add here an extra delete all files in folder command
            if (rank ==0){
                deleteall(newfolder);
            }
            ctx->barrier();
            runtime = 0;

            lmp=firstgo(ctx,structurefiles, Outputfolder, runid, leftbound, rightbound, runtype, direction,CUSTOMOP,op.get(),1,nsteps,trail, rootdirectory, calculator, strides,deltaT,constanttopology);
        }
        if (rank ==0){
            CleanUp(newfolder, runtype, direction, leftbound, rightbound, weshouldstop,nsteps, strides,deltaT,calculator, usePLMD,0);
            auto result = Runtimechecker(newfolder);
            runtime   = result.first/deltaT;
            propername = result.second;
        }//this is a rather long one make sure all things here are fine, must have runtime and propername transfer and have return statement if colvar doesn't develop
        ctx->barrier();
        ctx->bcast_longlong(runtime, 0);
        ctx->bcast_string(propername, 0);
        ctx->bcast_int(weshouldstop, 0);
        DestroyLammps(lmp);
        if (runtime<nsteps){
            TerminationStall(weshouldstop); //so another way we can terminate without plumed intervention causing weshouldstop to trigger termination stall is if we hit nsteps, this is a failsafe to allow the code to exit
        }
    }
}

template <typename T>// creating a prompter so as to not repeat pattern so often, have to make type fluid
T ask(const std:: string& prompt, std::function<bool(const T&)> validator = [] (const T&){return true;}){
    T value;
    while (true){
        std::cout<< prompt<<": ";
        std::cin >>value;
        if (!std::cin) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid"<<std::endl;
            continue;
        }

        if (!validator(value)) {
            std::cout << "Invalid"<<std::endl;
            continue;
        }

        return value;
    }
}


int isinteger(const std::string& s){
    if (s.empty()) return 0;
    for (char c:s){
        if (!std::isdigit(static_cast<unsigned char>(c))) return 0;
    }
    return 1;
}

int checkstructureinputs(std::string Delimeter, std::string structurefiles, int constanttopology){
    if (!std::filesystem::exists(structurefiles)){
        std::cout<<"Structure File Directory Doesn't Exist"<<std::endl;
        return 0;
    }
    if (!std::filesystem::is_directory(structurefiles)){
        std::cout<<"Structure File Directory is not a directory (is probably a file)";
        return 0;
    }
    int foundAtleastOne = 0;
    for (const auto& entry : std::filesystem::directory_iterator(structurefiles)){
        if (!entry.is_regular_file()){
            continue;
        }
        std::filesystem::path p = entry.path();
        if (p.extension().string() == Delimeter){
            foundAtleastOne = 1;
            std::string stem = p.stem().string();
            if (!isinteger(stem)){
                std::cout<<"Error: File "<<p.filename()<<" does not have an integer name (naming convention must be int.gro)"<<std::endl;
                return 0;
            }
            if (Delimeter == ".gro" and constanttopology ==0){
                std::filesystem::path ndx = std::filesystem::path(structurefiles) / (stem + ".ndx");
                std::filesystem::path top = std::filesystem::path(structurefiles) / (stem + ".top");
                if (!std::filesystem::exists(ndx)){
                    std::cout<<"Error: Missing "<<stem<<".ndx for "<<stem<<".gro"<<std::endl;
                    return 0;
                }
                if (!std::filesystem::exists(top)){
                    std::cout<<"Error: Missing "<<stem<<".top for "<<stem<<".gro"<<std::endl;
                    return 0;
                }
                if (!isinteger(ndx.stem().string())){
                    std::cout<<"Error: "<<ndx.filename()<<" does not have an integer name (naming convention must be int.ndx)"<<std::endl;
                    return 0;
                }
                if (!isinteger(top.stem().string())){
                    std::cout<<"Error: "<<top.filename()<<" does not have an integer name (naming convention must be int.top)"<<std::endl;
                    return 0;
                }
            }
        }
    }
    if (!foundAtleastOne){
        std::cout<<"Error: No files with delimeter "<<Delimeter<<" found in "<<structurefiles<<"\n";
        return 0;
    }
    return 1;
}
std::string safe(const std::string &s) {
    return s.empty() ? "\"\"" : s; //for empty string cases being printed into bash
}

//this function makes the bash scripts
void writebash(std::string CPartition, std::string Cruntime,std::string email, std::string RPartition, std::string Rruntime, int cores, int MDengine, std::string CUSTOMOP, std::string Outputfolder, int runtype, int stride, std::string Delimeter, std::string structurefilesloc, std::string rootdirectory, long long nsteps, int batchsize, int perfoldertrajectories, double leftbound, double rightbound, int maxcrossings,std::string gromacslocation, std::string plumedlocation, std::string OPENMPI, std::string FFTW,double deltaT, int constanttopology){
    std::string commandline ="srun -n 1 ./AdSamp " +std::to_string(runtype) + " "+ safe(Outputfolder) + " "+ std::to_string(leftbound) + " "+ std::to_string(rightbound) +" "+ std::to_string(maxcrossings);
    //adjust to also pass the batchsize and runs per trajectory
    std::ofstream controller("AdSampControl.sh");
    controller<<"#!/bin/bash"<<std::endl;
    controller<<"#SBATCH -p "<<CPartition<<std::endl;
    controller<<"#SBATCH -n 1"<<std::endl;
    controller<<"#SBATCH -N 1"<<std::endl;
    controller<<"#SBATCH -t "<<Cruntime<<std::endl;
    controller<<"#SBATCH --signal=B:10@60"<<std::endl;
    controller<<"#SBATCH --requeue"<<std::endl;
    if (email != "0"){
        controller<<"#SBATCH --mail-user="<<email<<std::endl;
        controller<<"#SBATCH --mail-type=end"<<std::endl;
    }

    controller<<"if [ -f /etc/bashrc ]; then"<<std::endl;
    controller<<"        . /etc/bashrc"<<std::endl;
    controller << "fi" << std::endl;
    controller<< "module purge"<<std::endl;
    controller<<"module load "<<OPENMPI<<std::endl;
    controller<<"module load "<<FFTW<<std::endl;
    controller<<"export PATH=$PATH:"<<plumedlocation<<"/bin:"<<gromacslocation<<"/bin"<<std::endl;
    controller<<"export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"<<plumedlocation<<"/lib:"<<gromacslocation<<"/lib"<<std::endl;
    controller<<"trap \"echo -n 'TIMEOUT @ '; date; echo 'Resubmitting...'; scontrol requeue ${SLURM_JOB_ID}\" 10"<<std::endl;
    controller<<commandline<<std::endl;
    controller.close();


    //add gromacs load if gromacs involved ugh, will have to prompt the user for its location
    std::ofstream individual("AdSamp.sh");
    individual<<"#!/bin/bash"<<std::endl;
    individual<<"#SBATCH -p "<<RPartition<<std::endl;
    individual<<"#SBATCH -n "<< cores<<std::endl;
    individual<<"#SBATCH -N 1"<<std::endl;
    individual<<"#SBATCH -t "<<Rruntime<<std::endl;
    individual<<"#SBATCH -a 1-"<<batchsize<<std::endl;
    individual<<"#SBATCH --signal=B:10@60"<<std::endl;
    individual<<"#SBATCH --requeue"<<std::endl;
    individual<<"trap \"echo -n 'TIMEOUT @ '; date; echo 'Resubmitting...'; scontrol requeue ${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}\" 10"<<std::endl;
    individual<<"if [ -f /etc/bashrc ]; then"<<std::endl;
    individual<<"        . /etc/bashrc"<<std::endl;
    individual << "fi" << std::endl;
    individual<<"module purge" <<std::endl;
    individual<<"module load "<<OPENMPI<<std::endl;
    individual<<"module load "<<FFTW<<std::endl;
    individual<<"export PATH=$PATH:"<<plumedlocation<<"/bin:"<<gromacslocation<<"/bin"<<std::endl;
    individual<<"export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"<<plumedlocation<<"/lib:"<<gromacslocation<<"/lib"<<std::endl;
    if (MDengine ==0){
        individual << "srun --cpu-bind=none -n " << cores<< " ./AdSamp "<<std::to_string(MDengine)<<" "<<std::to_string(runtype)<<" "<< safe(Outputfolder)<<" "<<"${SLURM_ARRAY_TASK_ID} "<<std::to_string(leftbound)<<" "<<std::to_string(rightbound)<<" "<<safe(structurefilesloc)<<" "<<Delimeter<<" "<<std::to_string(nsteps)<<" "<<safe(CUSTOMOP)<<" "<<safe(rootdirectory)<<" "<<std::to_string(perfoldertrajectories)<<" "<<std::to_string(stride)<<" "<<std::to_string(deltaT)<<" "<<std::to_string(constanttopology)<< " &"<< std::endl;
    }
    else{
        individual << "./AdSamp "<<std::to_string(MDengine)<<" "<<std::to_string(runtype)<<" "<< safe(Outputfolder)<<" "<<"${SLURM_ARRAY_TASK_ID} "<<std::to_string(leftbound)<<" "<<std::to_string(rightbound)<<" "<<safe(structurefilesloc)<<" "<<Delimeter<<" "<<std::to_string(nsteps)<<" "<<safe(CUSTOMOP)<<" "<<safe(rootdirectory)<<" "<<std::to_string(perfoldertrajectories)<<" "<<std::to_string(stride)<<" "<<std::to_string(deltaT)<<" "<<std::to_string(constanttopology)<< " &"<< std::endl;
    }
    individual<<"wait"<<std::endl;
    individual.close();
}

int gronamechecker(std::string rootdirectory, std::string structurefiles,std::string Delimeter, std::string Delimeter2){ //check if every structure file is mapped to an original file via gronames
    std::filesystem::path gronamespath =std::filesystem::path(structurefiles) / "gronames.txt";
    if (!std::filesystem::exists(gronamespath)){
        std::cout<<"Error, no gronames file in structures directory "<<structurefiles<<std::endl;
        return 0;
    }
    std::ifstream file(gronamespath);
    if (!file){
        std::cout<<"Can't open groname file"<<std::endl;
        return 0;
    }
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(file, line)){
        lines.push_back(line); //assemble all the entries in gronames.txt
    }
    file.close();

    for (const auto& entry : std::filesystem::directory_iterator(structurefiles)){//check for all files with delimeter
        if (!entry.is_regular_file()){
            continue;
        }
        std::filesystem::path p = entry.path();
        if (p.extension().string() == Delimeter){
            std::string stem = p.stem().string();
            int matchedline = 0;
            for (const std::string& line : lines){
                std::stringstream ss(line);
                std::string token;
                std::vector<std::string> tokens;
                while (std::getline(ss, token, ',')){
                    while (!token.empty() && token.front() == ' ') token.erase(token.begin());
                    while (!token.empty() && token.back() == ' ') token.pop_back();
                    tokens.push_back(token);
                }
                if (tokens.empty()){
                    continue;
                }
                if (tokens.back() == stem){
                    matchedline = 1;
                    std::size_t pos = line.find("initialconfiguration_");
                    if (pos == std::string::npos){
                        std::cout<<"Error line matching "<<p.filename().string()<<" in gronames.txt does not contain initialconfiguration_: "<<line<<std::endl;
                        return 0;
                    }
                    pos += std::string("initialconfiguration_").size();
                    std::size_t endpos = line.find('_', pos);
                    if (endpos == std::string::npos){
                        std::cout<<"Error could not get integer after initialconfiguration_ in line: "<<line<<std::endl;
                        return 0;
                    }
                    std::string numberstr = line.substr(pos, endpos - pos);
                    if (!isinteger(numberstr)){
                        std::cout<<"Error value after initialconfiguration_ is not an integer in line: "<<line<<std::endl;
                        return 0;
                    }
                    std::string expectedfile = std::to_string(std::stoi(numberstr)) + Delimeter2;
                    std::filesystem::path expectedpath = std::filesystem::path(rootdirectory) / expectedfile;
                    if (!std::filesystem::exists(expectedpath)){
                        std::cout<<"Error expected .gro file "<<expectedpath<<" does not exist for "<<p.filename().string()<<std::endl;
                        return 0;
                    }
                    break;
                }   
            }
            if (!matchedline){
                std::cout<<"Error: "<<p.filename().string()<<" does not appear as the last comma-separated entry in any line of gronames.txt"<<std::endl;
                return 0;
            }
        }
    }
    return 1;

}
void processingbash(std::string OPENMPI, std::string FFTW, std::string gromacslocation, std::string plumedlocation, std::string partition, std::string runtime, std::string email, int runtype, std::string runfolder, std::string newfoldername, std::string prior ){
    std::string commandline ="srun -n 1 ./AdSamp " +safe(runfolder) + " "+ safe(newfoldername) + " "+ safe(prior) + " "+ std::to_string(runtype);
    //adjust to also pass the batchsize and runs per trajectory
    std::ofstream controller("PostProcessing.sh");
    controller<<"#!/bin/bash"<<std::endl;
    controller<<"#SBATCH -p "<<partition<<std::endl;
    controller<<"#SBATCH -n 1"<<std::endl;
    controller<<"#SBATCH -N 1"<<std::endl;
    controller<<"#SBATCH -t "<<runtime<<std::endl;
    if (email != "0"){
        controller<<"#SBATCH --mail-user="<<email<<std::endl;
        controller<<"#SBATCH --mail-type=end"<<std::endl;
    }

    controller<<"if [ -f /etc/bashrc ]; then"<<std::endl;
    controller<<"        . /etc/bashrc"<<std::endl;
    controller << "fi" << std::endl;
    controller<<"module purge"<<std::endl;
    controller<<"module load "<<OPENMPI<<std::endl;
    controller<<"module load "<<FFTW<<std::endl;
    controller<<"export PATH=$PATH:"<<plumedlocation<<"/bin:"<<gromacslocation<<"/bin"<<std::endl;
    controller<<"export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"<<plumedlocation<<"/lib:"<<gromacslocation<<"/lib"<<std::endl;

    controller<<commandline<<std::endl;
}

void exportfromcptbash(std::string OPENMPI, std::string FFTW, std::string gromacslocation, std::string plumedlocation, std::string partition, std::string runtime, std::string email, int MDEngine, std::string postprocessedfolder, std::string newfoldername, int constanttopology, std::string rootdirectory){
    //adjust to also pass the batchsize and runs per trajectory
    std::ofstream controller("CPTExport.sh");
    controller<<"#!/bin/bash"<<std::endl;
    controller<<"#SBATCH -p "<<partition<<std::endl;
    controller<<"#SBATCH -n 1"<<std::endl;
    controller<<"#SBATCH -N 1"<<std::endl;
    controller<<"#SBATCH -t "<<runtime<<std::endl;
    if (email != "0"){
        controller<<"#SBATCH --mail-user="<<email<<std::endl;
        controller<<"#SBATCH --mail-type=end"<<std::endl;
    }

    controller<<"if [ -f /etc/bashrc ]; then"<<std::endl;
    controller<<"        . /etc/bashrc"<<std::endl;
    controller << "fi" << std::endl;
    controller << "module purge" << std::endl;
    controller<<"module load "<<OPENMPI<<std::endl;
    controller<<"module load "<<FFTW<<std::endl;
    controller<<"export PATH=$PATH:"<<plumedlocation<<"/bin:"<<gromacslocation<<"/bin"<<std::endl;
    controller<<"export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"<<plumedlocation<<"/lib:"<<gromacslocation<<"/lib"<<std::endl;
    
    if (MDEngine ==0){
        controller << "srun -n 1 ./AdSamp CONVERT " +safe(postprocessedfolder) + " "+ safe(rootdirectory) + " "+ safe(newfoldername) + " "+ std::to_string(MDEngine) + " " + std::to_string(constanttopology)<<std::endl; //add in an extra input so not same length as controller
    }
    else{
        controller << "./AdSamp CONVERT " +safe(postprocessedfolder) + " "+ safe(rootdirectory) + " "+ safe(newfoldername) + " "+ std::to_string(MDEngine) + " " + std::to_string(constanttopology)<<std::endl; //add in an extra input so not same length as controller
    }
}

void Login(){ //this will take in user inputs, ensure files are in order, then create a bash script for controller launch
    //part one packages
    std::string OPENMPI  = ask<std::string>("OpenMPI package name");
    std::string FFTW = ask<std::string>("FFTW package name");
    int MDengine = ask<int>("MD Calculator (0: LAMMPS, 1: GROMACS)",[](int x){return x==0 || x==1;});
    std::string gromacslocation="";
    if (MDengine ==1){
        gromacslocation = ask<std::string>("Gromacs Installation Location (ex: /home/mal258/ManuelsGromacs)",[](const std::string& s){return std::filesystem::exists(s+"/bin") && std::filesystem::is_directory(s+"/bin");});//check for bin
    }
    std::string plumedlocation = ask<std::string>("Plumed Installation Location (ex: /home/mal258/ManuelsPlumed)",[](const std::string& s){return std::filesystem::exists(s+"/bin") && std::filesystem::is_directory(s+"/bin");});
    int isPost= ask<int>("Postprocessing (Yes: 1, No: 2)",[](int x){return x==1 || x==2;});
    if (isPost ==1){
        int isExport = ask<int>("Exporting from binary files? (after postprocessing) (Yes: 1, No: 2)",[](int x){return x==1 || x==2;});
        std::string partition = ask<std::string>("Partition");
        std::string runtime = ask<std::string>("Max Run Time (D-HH:MM:SS)",[](const std::string& s){std::regex r(R"((\d+-)?\d{1,2}:\d{2}:\d{2})");return std::regex_match(s, r);});
        std::string email= ask<std::string>("email (write 0 for no email)");
        if (isExport ==2){
            int runtype = ask<int>("Most Recent Run type (0: Basin, 1: Probabilistic, 2: Exploration)",[](int x){return x==0 || x==1 || x ==2;});
            std::string runfolder = ask<std::string>("Completed Run Directory",[](const std::string& s){return std::filesystem::exists(s) && std::filesystem::is_directory(s);});
            std::string newfoldername = ask<std::string>("Postprocessing folder name",[](const std::string& s){return !std::filesystem::exists(s);});;
            std::string prior="";
            if (runtype ==1){
                prior = ask<std::string>("Groname File",[](const std::string& s){return std::filesystem::exists(s);});
            }
            processingbash(OPENMPI, FFTW, gromacslocation, plumedlocation, partition, runtime, email, runtype, runfolder, newfoldername, prior);
            return;
        }else{//export from cpt requires processedfolder, rootdirectory, outputfolder, MDengine, and constant topology
            int constanttopology=0;
            if (MDengine ==1){
                constanttopology = ask<int>("Constant Topology? (one index and topology file in working directory applies to all structures), 0: no, 1: yes",[](const int x){if (x != 0 && x != 1) return false;if (x == 1){return std::filesystem::exists("topol.top") && std::filesystem::exists("index.ndx");}return true;});
            }
            int condition = 0;
            std::string postprocessedfolder;
            while (condition ==0){
                postprocessedfolder = ask<std::string>("Postprocessed Directory",[](const std::string& s){return std::filesystem::exists(s) && std::filesystem::is_directory(s);});
                condition = checkstructureinputs(".cpt", postprocessedfolder, constanttopology);
            }
            std::string newfoldername = ask<std::string>("Output folder name",[](const std::string& s){return !std::filesystem::exists(s);});;
            condition = 0;
            std::string rootdirectory;
            if (MDengine ==1){ //lammps doesn't need this matching
                while (condition==0){
                    rootdirectory = ask<std::string>("Root directory (where .gro files are)",[](const std::string& s){return std::filesystem::exists(s) && std::filesystem::is_directory(s);});
                    condition = checkstructureinputs(".gro", rootdirectory,constanttopology); //reusing the same function by changing delimeter
                    if (condition!=0){
                        condition = gronamechecker(rootdirectory, postprocessedfolder, ".cpt", ".gro"); //for several ffs steps postprocessing will create a gronames file that gives the original .gro file
                    }
                }
            }
            exportfromcptbash(OPENMPI, FFTW, gromacslocation, plumedlocation, partition, runtime, email, MDengine, postprocessedfolder, newfoldername, constanttopology, rootdirectory);
            return;
        }
    }
    //part 2 operating system
    std::string CPartition = ask<std::string>("Controller Partition (maximize runtime)");
    std::string Cruntime = ask<std::string>("Controller Max Run Time (D-HH:MM:SS)",[](const std::string& s){std::regex r(R"((\d+-)?\d{1,2}:\d{2}:\d{2})");return std::regex_match(s, r);});
    std::string email= ask<std::string>("email (write 0 for no email)");
    std::string RPartition = ask<std::string>("Runner Partition (maximize cores)");
    std::string Rruntime = ask<std::string>("Runner Max Run Time (D-HH:MM:SS)",[](const std::string& s){std::regex r(R"((\d+-)?\d{1,2}:\d{2}:\d{2})");return std::regex_match(s, r);});
    int cores = ask<int>("Runner Cores",[](int x){return x>0;});
    //run info
    double deltaT = ask<double>("MD Step Size (ps)",[](double x){return x>0;});    
    //preemptive mdp check
    if (!std::filesystem::exists("ffs.mdp")){
        std::cout<<"Missing ffs.mdp file in working directory"<<std::endl; //I realize this takes up a lot of lines and is kinda a convoluted part but better to fail early rather than cause an issue on a cluster
    }
    std::ifstream mdpfile("ffs.mdp");
    std::string line;
    int hasvel=0;
    int hastemp=0;
    int hasseed=0;
    int hasnsteps =0;
    int hasdelT=0;
    if (MDengine ==1){
        while (std::getline(mdpfile, line)){
            if (line.find("gen-vel")!=std::string::npos) {
                hasvel = 1;
            }
            if (line.find("gen-temp")!=std::string::npos){
                hastemp=1;
            }
            if (line.find("gen-seed")!= std::string::npos){
                hasseed =1;
            }
            if (line.find("nsteps")!= std::string::npos){
                hasnsteps = 1;
            }
            if (line.find("dt ") != std::string::npos){
                hasdelT = 1;
            }
        }
        if (hasvel ==0 || hastemp ==0 || hasseed ==0 || hasnsteps ==0 || hasdelT ==0){
            std::cout<<"Error, faulty ffs.mdp"<<std::endl;
            if (hasvel ==0){
                std::cout<<"Missing gen-vel line"<<std::endl;
            }if (hastemp ==0){
                std::cout<<"Missing gen-temp line"<<std::endl;
            }if (hasseed==0){
                std::cout<<"Missing gen-seed line"<<std::endl;
            }if (hasnsteps==0){
                std::cout<<"Missing nsteps"<<std::endl;
            }if (hasdelT==0){
                std::cout<<"Missing dt"<<std::endl;
            }
            return;
        }
    }else{
        while (std::getline(mdpfile, line)){
            if ((line.find("npt")!=std::string::npos && line.find("temp")!=std::string::npos)|| (line.find("nvt") !=std::string::npos && line.find("temp")!=std::string::npos)){
                hasvel =1;
            }
            if (line.find("timestep")!= std::string::npos){
                hasdelT =1;
            }
            if (line.find("velocity")!=std::string::npos){
                std::cout<<"remove velocity initialization line, code has internal method"<<std::endl;
                return;
            }
        }
        if (hasvel ==0 || hasdelT ==0){
            std::cout<<"Error, faulty ffs.mdp"<<std::endl;
            if (hasvel==1){
                std::cout<<"Missing temperature (npt/nvt) line"<<std::endl;
            }
            if (hasdelT==0){
                std::cout<<"Missing timestep line"<<std::endl;
            }
            return;
        }
    }
    std::string CUSTOMOP;
    int usePLMD=1;
    if (MDengine ==0){
        usePLMD = ask<int>("Using PLUMED (1: yes, 0: no)",[](int x){return x==1 || x==0;});
        if (usePLMD ==0){
            CUSTOMOP = ask<std::string>("Custom OP Filename");
        }
    }
    if (usePLMD ==1){
        if (!std::filesystem::exists("plumed.dat")){
            std::cout<<"Missing plumed.dat file in working directory"<<std::endl;
            return;
        }
        std::ifstream plumedfile("plumed.dat");
        std::string line;
        int hascommitor=0;
        int hasprint=0;
        int hasflush=0;
        int tick =0;
        while (std::getline(plumedfile, line)){
            if ((line.find("ARG=") != std::string::npos && tick ==1) || (line.find("STRIDE=") != std::string::npos && tick ==2) || (line.find("BASIN_LL1") != std::string::npos && tick ==3) || (line.find("BASIN_UL1") != std::string::npos && tick ==4) || (line.find("BASIN_LL2") != std::string::npos && tick ==5) || (line.find("BASIN_UL2") != std::string::npos && tick ==6) || (line.find("... COMMITTOR") != std::string::npos && tick ==7) ){
                tick++;
            }
            if (line.find("COMMITTOR") != std::string::npos && line.find("... COMMITTOR") == std::string::npos){
                hascommitor = 1;
                tick +=1;
            }
            if (line.find("PRINT")!=std::string::npos){
                hasprint=1;
            }
            if (line.find("FLUSH")!= std::string::npos){
                hasflush =1;
            }
            if (hascommitor ==1 && tick ==0){
                tick +=1;
            }
        }
        if (hasprint ==0 || hasflush ==0 || hascommitor ==0 || tick!=8){
            if (hasprint==0){
                std::cout<<"Missing Print Line"<<std::endl;
            }
            if (hasflush==0){
                std::cout<<"Missing Flush Line"<<std::endl;
            }
            if (hascommitor==0){
                std::cout<<"Missing Commitor Line"<<std::endl;
            }
            if (tick !=8){
                std::cout<<"Commitor Section Formatting Is Incorrect (refer to manual)"<<std::endl;
            }
            std::cout<<"Error, faulty plumed.dat"<<std::endl;
            return;
        }
    }
    std::string Outputfolder = ask<std::string>("Output folder",[](const std::string & s){static const std::regex r(R"(^[A-Za-z0-9._-]+$)"); return std::regex_match(s, r);});
    int runtype = ask<int>("Run type (0: Basin, 1: Probabilistic, 2: Exploration)",[](int x){return x==0 || x==1 || x ==2;});
    int stride = ask<int>("OP Stride",[](int x){return x>0;});
    std::string Delimeter;
    std::string structurefilesloc;
    int condition = 0;
    int tr=0;
    int constanttopology=0;
    while (condition ==0){ //condition is whether there are that structure file type in the folder also if runtype == .gro make sure .top and .ndx in there too
        tr=0;
        Delimeter = ask<std::string>("Structure File Delimeter (Including .)",[](const std::string& s){return s.find('.')!=std::string::npos;});
        structurefilesloc = ask<std::string>("Structure File Directory",[](const std::string& s){return std::filesystem::exists(s) && std::filesystem::is_directory(s);});
        if (MDengine==1){
            constanttopology = ask<int>("Constant Topology? (one index and topology file in working directory applies to all structures), 0: no, 1: yes",[](const int x){if (x != 0 && x != 1) return false;if (x == 1){return std::filesystem::exists("topol.top") && std::filesystem::exists("index.ndx");}return true;});
        }
        if (MDengine==1){
            if (Delimeter != ".gro" && Delimeter != ".cpt"){
                std::cout<<" For Gromacs Operation delimeter must be either .gro or .cpt "<<std::endl;
                condition =0;
                tr=1;
            }
        }
        if (tr!=1){
            condition = checkstructureinputs(Delimeter, structurefilesloc, constanttopology);
        }
    }
    std::string rootdirectory;
    if (Delimeter==".cpt" && MDengine ==1){
        condition = 0;
        while (condition==0){
            rootdirectory = ask<std::string>("Root directory (where .gro files are)",[](const std::string& s){return std::filesystem::exists(s) && std::filesystem::is_directory(s);});
            condition = checkstructureinputs(".gro", rootdirectory,constanttopology); //reusing the same function by changing delimeter
            if (condition!=0){
               condition = gronamechecker(rootdirectory, structurefilesloc, Delimeter, ".gro"); //for several ffs steps postprocessing will create a gronames file that gives the original .gro file
            }
        }
    }
    long long nsteps=1;
    int batchsize;
    int perfoldertrajectories=1;
    int maxcrossings = std::numeric_limits<int>::max();
    if (runtype!=1){
        nsteps = ask<int>("Desired Number of Steps",[](int x){return x>0;});
        batchsize = ask<int>("Batch Size (for 1 per file return 0)",[](int x){return x>-1;});
        if (batchsize ==0){
            std::vector<std::string> collectionofinputs = listdirEnding(structurefilesloc, Delimeter);
            batchsize = collectionofinputs.size();
        }
    }else{
        batchsize = ask<int>("Batch Size",[](int x){return x>-1;});
        maxcrossings = ask<int>("Max Crossings Before Interruption",[](int x){return x>-1;});
    }

    double leftbound=0; double rightbound=0;
    if (runtype< 2){
        while (leftbound==rightbound){
            leftbound = ask<double>("lambda_A");
            rightbound = ask<double>("Goal Lambda");
            if (leftbound==rightbound){
                std::cout<<"error, leftbound cant equal rightbound"<<std::endl;
            }
        }
    }
    if (runtype ==1){
        perfoldertrajectories = ask<int>("Trajectories per folder",[](int x){return x>0;});
    }
    writebash(CPartition, Cruntime, email, RPartition, Rruntime, cores, MDengine, CUSTOMOP, Outputfolder, runtype, stride, Delimeter, structurefilesloc, rootdirectory, nsteps, batchsize, perfoldertrajectories, leftbound, rightbound, maxcrossings,gromacslocation,plumedlocation,OPENMPI, FFTW,deltaT, constanttopology);
}
std::string getSlurmJobID() {
    const char* id = std::getenv("SLURM_JOB_ID");
    return id ? std::string(id) : "";
}

bool isJobRunning(int jobid) {
    std::string cmd = "squeue -h -r -j " + std::to_string(jobid) + " -o \"%i %R\"";//this tests if a job is currently running
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return false;
    char buffer[256];
    std::string target = std::to_string(jobid);
    int running = 0;
    while (fgets(buffer, sizeof(buffer), pipe)) {
        std::string line(buffer);
        if (!line.empty() && line.back() == '\n')
            line.pop_back();
        if (line == target || line.rfind(target + "_", 0) == 0) {
            if (line.find("launch failed requeued held") != std::string::npos) {
                std::cout<<"scontrol release " + line.substr(0, line.find(' '))<<std::endl;
                std::system(("scontrol release " + line.substr(0, line.find(' '))).c_str()); //Grace partition has a job launch held problem gonna resolve it by releasing it
            }else{
                running =1;
            }
        }
    }
    pclose(pipe);
    return running==1;  //true if one is properly running false otherwise
}
int submitjob(std::string scriptname){
    std::string cmd = "sbatch " + scriptname;
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) {
        std::cout<<"Job Submit Failure"<<std::endl;
        return 0;
    }
    char buffer[256];
    std::string output;
    while (fgets(buffer, sizeof(buffer), pipe)) {
        output += buffer;
    }
    pclose(pipe);
    //it might come to be the case that someone not on my tested clusters (grace, bouchet, expanse) will have different output when submitting a job
    //in manual flag this set of lines for edits
    std::regex re("Submitted batch job ([0-9]+)");
    std::smatch match;
    if (std::regex_search(output, match, re)) {
        return std::stoi(match[1]);
    }
    return 0;
}

int readbasinpathline(const std::string& line){
    size_t start = 0;
    for ( int i = 0; i<3; ++i){
        start = line.find(',', start);
        if (start ==std::string::npos){//here we check for spots to split by ',' and return 0 crossings if there aren't four 
            return 0;
        }
        start +=1; //keep going past the comma
    }
    size_t end = line.find(',', start);
    std::string token;
    if (end == std::string::npos){
        token = line.substr(start);
    }else{
        token = line.substr(start, end-start);
    }
    return std::stoi(token);
}

std::pair<std::string, std::string> DateandTime(){
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    std::tm* tm_ptr = std::localtime(&now_time);

    std::string str1 = std::to_string(tm_ptr->tm_mon + 1) + "/" +
                       std::to_string(tm_ptr->tm_mday) + "/" +
                       std::to_string(tm_ptr->tm_year + 1900);

    std::ostringstream time_stream;
    time_stream << tm_ptr->tm_hour << ":";
    if (tm_ptr->tm_min < 10) time_stream << "0";
    time_stream << tm_ptr->tm_min << ":";
    if (tm_ptr->tm_sec < 10) time_stream << "0";
    time_stream << tm_ptr->tm_sec;

    std::string str2 = time_stream.str();

    return {str1, str2};
}

int checkcrossings(int simtype, std::string outputfolder, double leftbound, double rightbound){
    int direction =rightbound>leftbound?1:-1; //pos 1 for increasing, neg one for decreasing
    if (!std::filesystem::exists("RunInfo.txt")){ //write headers
        std::ofstream control("RunInfo.txt");
        if (simtype ==0){ //basin
            control<<"Date, Time, Started Simulations, Crossings, Average Crossings Per Initial Configuration"<<std::endl;
        }else if(simtype ==1){ //probabilistic
            control<< "Date, Time, Started Simulations, Successes, Failures, Perc Success, ";
            if (direction==1){
                control<<"Max λ, ";
            }else{
                control<<"Min λ, ";
            }
            control<<"Ave Success Length, Average Failure Length"<<std::endl;
        }else if(simtype ==2){ //exploration
            control<<"Date, Time, Started Simulations, Completed Simulations"<<std::endl;
        }
        control.close();
    }
    //step 1: pipe in as safe of a way as possible all the paths within all the subfolders
    if (!std::filesystem::exists(outputfolder) || !std::filesystem::is_directory(outputfolder)){
        return 0;
    }
    std::vector<std::string> pathfiles;
    for (const auto&entry : std::filesystem::directory_iterator(outputfolder)){
        if (entry.is_directory()){
            std::filesystem::path p =entry.path() /"path.txt";//this gives the location of my path file
            if (std::filesystem::exists(p)){
                pathfiles.push_back(p.string());
            }
        }
    }
    if (pathfiles.empty()){
        return 0;
    }
    //gonna follow the pattern I used in the python FFS script, cat seems to be a safe way to treat the paths
    std::string cmd = "cat ";
    for (std::string& f : pathfiles) {
        cmd += "\"" + f + "\" ";
    }
    cmd += "> globalpath.txt";
    std::system(cmd.c_str());
    //step 2: now read these line by line to get the run paremeters

    std::ifstream path("globalpath.txt");
    std::vector<std::string> pathentries;
    int total=0;
    int successes =0;
    int failures = 0;
    std::string line;
    double stime=0;
    double ftime = 0;
    std::string propername;
    double maxval = -1* direction * std::numeric_limits<double>::max(); // big number negative for plus 1 pos otherwise
    while (std::getline(path, line)){
        total +=1;
        if (simtype ==0){
            successes += readbasinpathline(line);
        }if (simtype ==1){
            if (line.find("success")!= std::string::npos || line.find("failure")!= std::string::npos){
                pathentries = pathlinereader(line);
                propername = pathentries[2];
                std::string clave = "subfolder_"; // key 
                size_t pos = propername.find(clave);
                pos += clave.length();
                size_t end = pos;
                while (end<propername.size() && std::isdigit(propername[end])){
                    ++end;
                }
                std::string Colvar = outputfolder + "/" + std::to_string(std::stoi(propername.substr(pos, end -pos))) +"/C" + propername;
                auto result = lastorderparameter(Colvar);
                double time = result.first;
                double c= result.second;
                if (line.find("success")!= std::string::npos){
                    successes +=1;
                    if (c*direction>maxval*direction){
                        maxval = c;
                    }
                    stime += time;
                }else{
                    failures +=1; //don't have to go through the find fail rigamaroll cus of process of elimination
                    ftime +=time;
                }
            }
        }if (simtype ==2){//Exploration
            if (line.find("complete") != std::string::npos){
                successes +=1;
            }
        }
    }
    path.close();

    //now append to the controller file
    auto [date, time] = DateandTime();
    std::ofstream control("RunInfo.txt", std::ios::app);
    control<<date<<", "<<time<<", ";
    if (simtype ==0){
        double ave=0;
        if (successes !=0){
            ave = static_cast<double>(successes)/total; //alright it's pretty cool that c++ let's me control parallelism and has pointers and scopes but its terrible that this is how you get a double from integer division
        }
        control<<total<<", "<<successes<< ", "<< ave<<std::endl;
    }else if (simtype ==1){
        double percentsucc = 0;
        double avestime = 0;
        double aveftime = 0;
        if (successes!=0 ){
            avestime =static_cast<double>(stime)/successes;
            percentsucc = 1;
            if (failures!=0){
                percentsucc = 100.0 *successes/ (successes + failures); // the 100 is the double here
            }
        }else{
            maxval = 0;
        }
        if (failures !=0){
            aveftime = static_cast<double>(ftime)/failures;
        }
        control<<total<<", "<<successes<<", "<<failures<<", "<<percentsucc<<", "<< maxval<<", "<<avestime<<", "<<aveftime<<std::endl;
    }else{
        control<<total<<", "<<successes<<std::endl;
    }
    control.close();
    return successes;
}
void controller(int simtype, std::string outputfolder, double leftbound, double rightbound, int maxcrossings){
    if (!std::filesystem::exists(outputfolder)){
        std::filesystem::create_directories(outputfolder);
    }
    //controller steps 1. see if there's another active controller, if so check slurm tag to ensure this one is more recent, if so cancel it
    int prevjob=0;
    std::ifstream controljobs("controljobs.txt");
    if (controljobs.is_open()){ 
        if (!(controljobs >> prevjob)) {
            prevjob =0;
        }
    }
    controljobs.close();
    int myjob = std::stoi(getSlurmJobID());
    if (prevjob!=0){
        if (myjob>prevjob){
            std::cout<<"scancel " + std::to_string(prevjob)<<std::endl;
            std::string cmd = "scancel " + std::to_string(prevjob);
            std::system(cmd.c_str()); 
        }
    }
    std::ofstream out("controljobs.txt",std::ios::trunc);
    out<<myjob<<std::endl;
    out.close();
    //step 2. See if there are subjobs running, if there are, skip to controlling
    int subjobs=0;
    std::ifstream arrayjobs("arrayjobs.txt");
    if (arrayjobs.is_open()){
        if (!(arrayjobs>>subjobs)){
            subjobs = 0;
        }
    }
    arrayjobs.close();
    //step 3. if there aren't start jobs and record their job id
    if (subjobs ==0 || !isJobRunning(subjobs)){
        std::cout<<subjobs<<isJobRunning(subjobs)<<std::endl;
        subjobs = submitjob("AdSamp.sh");
        std::ofstream arrayjobs("arrayjobs.txt",std::ios::trunc);
        arrayjobs<<subjobs<<std::endl;
        arrayjobs.close();
    }
    //step 4. update controller outputs regularly for crossings stop once ncrossings>maxcrossings or sub jobs stop running 
    int ncrossings = 0;
    while (isJobRunning(subjobs) && ncrossings<maxcrossings){
        ncrossings = checkcrossings(simtype, outputfolder, leftbound, rightbound);
        std::this_thread::sleep_for(std::chrono::seconds(60)); //so that we're not checking so so so often
    }
    ncrossings = checkcrossings(simtype, outputfolder, leftbound, rightbound);
    std::cout<<isJobRunning(subjobs)<<std::endl;
    std::cout<<ncrossings<<" "<<maxcrossings<<std::endl;
    //std::string cmd = "scancel " + std::to_string(subjobs);
    //std::system(cmd.c_str()); 
    
}

double extractTime(const std::string& line) {
    const std::string token = "Time: ";
    size_t pos = line.find(token);
    return std::stod(line.substr(pos + token.size())); //extract from BasinCrossingRecorder.txt
}

void PostProcessing(std::string RUNFOLDER, std::string NEWFOLDER, std::string prior, int simtype){ 
    if (!std::filesystem::is_directory(NEWFOLDER)){
        std::filesystem::create_directories(NEWFOLDER);
    }
    if (std::filesystem::is_regular_file("RunInfo.txt")){
        std::filesystem::rename("RunInfo.txt", NEWFOLDER + "/RunInfo.txt");
    }
    if (simtype ==0){ //basin case
        int crossing=0;
        //in a basin run we should: 1. copy over the controller, 2. copy path lines , 3. copy all colvar files, 4. copy all checkpoint files and name them per global crossings, 5. give a mapping of crossing to file to time
        for (const auto& entry : std::filesystem::directory_iterator(RUNFOLDER)){
            std::string subfolder = entry.path().filename().string();
            //step 1, open path and find out the name
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/path.txt")){
                std::ifstream path(RUNFOLDER + "/" + subfolder + "/path.txt");
                std::string line;
                std::getline(path, line); //there should only be one line for basin
                std::ofstream allpath(NEWFOLDER + "/path.txt",std::ios::app);
                allpath<<line<<std::endl;
                allpath.close();
                std::vector<std::string> pathentries = pathlinereader(line);
                std::string propername = pathentries[2];
                if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/C" + propername )){
                    std::filesystem::copy_file(RUNFOLDER + "/" + subfolder + "/C" + propername, NEWFOLDER + "/C" + propername);
                }
                std::string ic = "initialconfiguration_";
                size_t start =propername.find(ic); //need to extract ic from propername
                start += ic.length();
                size_t end = propername.find("_", start);
                int chosenstructure = std::stoi(propername.substr(start,end-start));
                int ncrossings=0;
                if (pathentries.size()>3){
                    ncrossings = std::stoi(pathentries[3]);
                }
                int thisfolderscrossing=0;
                if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/BasinCrossingRecorder.txt")){
                    //want to create an output file that says original structure, crossing time, crossing number 
                    std::ifstream BasinCrossings(RUNFOLDER + "/" + subfolder + "/BasinCrossingRecorder.txt");
                    std::string Bline;
                    while (std::getline(BasinCrossings, Bline)){
                        crossing++;
                        thisfolderscrossing++;
                        double time = extractTime(Bline);
                        std::ofstream AllCrossings(NEWFOLDER + "/AllCrossings.txt", std::ios::app);
                        AllCrossings<< chosenstructure<< " , "<< propername<< " , "<< time<< " , "<< crossing<< std::endl;
                        AllCrossings.close();
                        if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/" + std::to_string(thisfolderscrossing) + ".cpt")){
                            std::filesystem::copy_file(RUNFOLDER + "/" + subfolder + "/" + std::to_string(thisfolderscrossing) + ".cpt", NEWFOLDER + "/" + std::to_string(crossing) + ".cpt");
                        }
                        std::ofstream newgronames(NEWFOLDER + "/gronames.txt", std::ios::app);
                        newgronames << propername<< " , "<<crossing<<std::endl;
                        newgronames.close();
                    }
                }
            }
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/outputfiles.txt")){
                std::filesystem::copy_file(RUNFOLDER +"/" + subfolder + "/outputfiles.txt",NEWFOLDER + "/folder_" + subfolder + "_outputfiles.txt");
            }
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/seedvalues.txt")){
                std::filesystem::copy_file(RUNFOLDER +"/" + subfolder + "/seedvalues.txt",NEWFOLDER + "/folder_" + subfolder + "_seedvalues.txt");
            }
        
        }
    }if (simtype ==1){ //probabilistic case
        std::unordered_map<int , std::string> PrevGro;
        if (!prior.empty()){
            std::ifstream gronames(prior);
            std::string line;
            while (std::getline(gronames, line)){
                size_t pos = line.rfind(',');
                PrevGro[std::stoi(line.substr(pos+1))] = line;
            }
        }
        int successcount = 0;
        for (const auto& entry : std::filesystem::directory_iterator(RUNFOLDER)){
            std::string subfolder = entry.path().filename().string();
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/path.txt")){
                std::ifstream path(RUNFOLDER + "/" + subfolder + "/path.txt");
                std::string line;
                while (std::getline(path, line)){
                    std::ofstream allpath(NEWFOLDER + "/path.txt", std::ios::app);
                    allpath<<line<<std::endl;
                    allpath.close();
                    std::vector<std::string> pathentries = pathlinereader(line);
                    std::string propername = pathentries[2];
                    if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/C" + propername )){
                        std::filesystem::copy_file(RUNFOLDER + "/" + subfolder + "/C" + propername, NEWFOLDER + "/C" + propername);
                    }
                    std::string ic = "initialconfiguration_";
                    size_t start =propername.find(ic); //need to extract ic from propername
                    start += ic.length();
                    size_t end = propername.find("_", start);
                    int chosenstructure = std::stoi(propername.substr(start,end-start));
                    if (line.find("success")!= std::string::npos){//success case
                        successcount++;
                        std::ofstream newgronames(NEWFOLDER + "/gronames.txt", std::ios::app);
                        if (!prior.empty()){
                            std::string prevline = PrevGro[chosenstructure];
                            newgronames<< prevline<< " , ";
                        }
                        newgronames << propername<< " , "<<successcount<<std::endl; // this way we have the proper output for both cases of first gronames and nto
                        //if we have a success we must copy the cpt file, copy the colvar file, append to successes.txt write last op to orderparameters.txt, add to gronames.txt
                        if (std::filesystem::is_regular_file(RUNFOLDER +"/"+ subfolder + "/" + propername + ".cpt")){
                            std::filesystem::copy_file(RUNFOLDER + "/" + subfolder + "/" + propername + ".cpt", NEWFOLDER + "/" + std::to_string(successcount) + ".cpt");
                        }
                        if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/C" + propername)){
                            auto result = lastorderparameter(RUNFOLDER + "/" + subfolder + "/C" + propername);
                            double time = result.first;
                            double c= result.second;
                            std::ofstream successestxt(NEWFOLDER + "/successes.txt", std::ios::app);
                            successestxt<<propername <<" , "<< c<<std::endl;
                            successestxt.close();
                            std::ofstream orderparameterstxt(NEWFOLDER + "/orderparameters.txt", std::ios::app);
                            orderparameterstxt<<propername<<" , "<<c<<std::endl;
                            orderparameterstxt.close();
                        }else{
                            std::ofstream successestxt(NEWFOLDER + "/successes.txt", std::ios::app);
                            successestxt<<propername <<" , NA"<<std::endl;;
                            successestxt.close();
                            std::ofstream orderparameterstxt(NEWFOLDER + "/orderparameters.txt", std::ios::app);
                            orderparameterstxt<<propername<<" , NA"<<std::endl;
                            orderparameterstxt.close();
                        }

                    }else if (line.find("failure")!=std::string::npos){
                        if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/C" + propername)){
                            auto result = lastorderparameter(RUNFOLDER + "/" + subfolder + "/C" + propername);
                            double time = result.first;
                            double c= result.second;
                            std::ofstream failurestxt(NEWFOLDER + "/failures.txt", std::ios::app);
                            failurestxt<<propername <<" , "<< c<<std::endl;
                            failurestxt.close();
                            std::ofstream orderparameterstxt(NEWFOLDER + "/orderparameters.txt",std::ios::app);
                            orderparameterstxt<<propername<<" , "<<c<<std::endl;
                            orderparameterstxt.close();
                        }else{
                            std::ofstream failurestxt(NEWFOLDER + "/failures.txt", std::ios::app);
                            failurestxt<< propername <<" , NA"<<std::endl;
                            failurestxt.close();
                            std::ofstream orderparameterstxt(NEWFOLDER + "/orderparameters.txt", std::ios::app);
                            orderparameterstxt<<propername<<" , NA"<<std::endl;
                            orderparameterstxt.close();
                        }
                    }else{ //incompletes
                        if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/C" + propername)){
                            std::ifstream checkFile(RUNFOLDER + "/"+ subfolder + "/C" + propername);
                            int linecount = 0;
                            std::string checkline;
                            while (std::getline(checkFile, checkline)) linecount++;
                            checkFile.close();
                            if (linecount>0){
                                auto result = lastorderparameter(RUNFOLDER + "/" + subfolder + "/C" + propername);
                                double time = result.first;
                                double c= result.second;
                                std::ofstream incompletestxt(NEWFOLDER + "/incompletes.txt", std::ios::app);
                                incompletestxt<<propername <<" , "<< c<<std::endl;
                                incompletestxt.close();
                                std::ofstream orderparameterstxt(NEWFOLDER + "/orderparameters.txt", std::ios::app);
                                orderparameterstxt<<propername<<" , "<<c<<std::endl;
                                orderparameterstxt.close();
                            }else{
                                std::ofstream incompletestxt(NEWFOLDER + "/incompletes.txt", std::ios::app);
                                incompletestxt<<propername <<" , NA"<<std::endl;;
                                incompletestxt.close();
                                std::ofstream orderparameterstxt(NEWFOLDER + "/orderparameters.txt",std::ios::app);
                                orderparameterstxt<<propername<<" , NA"<<std::endl;
                                orderparameterstxt.close();
                            }
                        }else{
                            std::ofstream incompletestxt(NEWFOLDER + "/incompletes.txt", std::ios::app);
                            incompletestxt<<propername <<" , NA"<<std::endl;;
                            incompletestxt.close();
                            std::ofstream orderparameterstxt(NEWFOLDER + "/orderparameters.txt",std::ios::app);
                            orderparameterstxt<<propername<<" , NA"<<std::endl;
                            orderparameterstxt.close();
                        }
                    }
                }
            }//remaining: copy iterationtracker.txt, outputfiles.txt, seedvalues.txt
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/iterationtracker.txt")){
                std::filesystem::copy_file(RUNFOLDER +"/" + subfolder + "/iterationtracker.txt",NEWFOLDER + "/folder_" + subfolder + "_iterationtracker.txt");
            }
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/outputfiles.txt")){
                std::filesystem::copy_file(RUNFOLDER +"/" + subfolder + "/outputfiles.txt",NEWFOLDER + "/folder_" + subfolder + "_outputfiles.txt");
            }
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/seedvalues.txt")){
                std::filesystem::copy_file(RUNFOLDER +"/" + subfolder + "/seedvalues.txt",NEWFOLDER + "/folder_" + subfolder + "_seedvalues.txt");
            }
        }

    }if (simtype ==2){ //exploration case
        //steps for exploration: copy all colvar files, check for complete, if complete copy cpt file
        for (const auto& entry : std::filesystem::directory_iterator(RUNFOLDER)){
            std::string subfolder = entry.path().filename().string();
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/path.txt")){
                std::ifstream path(RUNFOLDER + "/" + subfolder + "/path.txt");
                std::string line;
                std::getline(path, line); //there should only be one line for basin
                std::ofstream allpath(NEWFOLDER + "/path.txt",std::ios::app);
                allpath<<line<<std::endl;
                allpath.close();
                std::vector<std::string> pathentries = pathlinereader(line);
                std::string propername = pathentries[2];
                if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/C" + propername )){
                    std::filesystem::copy_file(RUNFOLDER + "/" + subfolder + "/C" + propername, NEWFOLDER + "/C" + propername);
                }
                std::string ic = "initialconfiguration_";
                size_t start =propername.find(ic); //need to extract ic from propername
                start += ic.length();
                size_t end = propername.find("_", start);
                int chosenstructure = std::stoi(propername.substr(start,end-start));
                
                if (line.find("complete") != std::string::npos){
                    //completed run case copy its checkpoint file and add its name to completedruns.txt
                    if (std::filesystem::is_regular_file(RUNFOLDER +"/"+ subfolder + "/" + propername + ".cpt")){
                        std::filesystem::copy_file(RUNFOLDER + "/" + subfolder + "/" + propername + ".cpt", NEWFOLDER + "/" + std::to_string(chosenstructure) + ".cpt");
                    }
                    std::ofstream completedruns(NEWFOLDER +"/completedruns.txt", std::ios::app);
                    completedruns<<propername<<std::endl;
                    completedruns.close();
                    std::ofstream newgronames(NEWFOLDER + "/gronames.txt", std::ios::app);
                    newgronames << propername<< " , "<<chosenstructure<<std::endl;
                    newgronames.close();
                }else{
                    std::ofstream incompleteruns(NEWFOLDER +"/incompleteruns.txt", std::ios::app);
                    incompleteruns<<propername<<std::endl;
                    incompleteruns.close();
                }
            }
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/outputfiles.txt")){
                std::filesystem::copy_file(RUNFOLDER +"/" + subfolder + "/outputfiles.txt",NEWFOLDER + "/folder_" + subfolder + "_outputfiles.txt");
            }
            if (std::filesystem::is_regular_file(RUNFOLDER + "/" + subfolder + "/seedvalues.txt")){
                std::filesystem::copy_file(RUNFOLDER +"/" + subfolder + "/seedvalues.txt",NEWFOLDER + "/folder_" + subfolder + "_seedvalues.txt");
            }
        }

    }
}
void generateimdp(){
    std::ifstream mdp("ffs.mdp");
    std::vector<std::string> lines;
    std::string line;
    int nstxoutline = -1;
    int nstvoutline = -1;
    int nstfoutline = -1;
    int nstepsline = -1; 
    int lineIndex = 0;
    while (std::getline(mdp,line)){
        auto pos = line.find("nstxout");
        if (pos != std::string::npos) {
            size_t after = pos + 7; 
            if (after >= line.size() || (!std::isalnum(line[after]) && line[after] != '-')) {
                nstxoutline = lineIndex;
            }
        }
        if (line.find("nstvout") !=std::string::npos){
            nstvoutline = lineIndex;
        }
        if (line.find("nstfout")!= std::string::npos){
            nstfoutline = lineIndex;
        }
        if (line.find("nsteps ") != std::string::npos){
            nstepsline = lineIndex;
        }
        lines.push_back(line);
        lineIndex++;
    } 
    mdp.close();
    if (nstxoutline!=-1 && nstvoutline!=-1 && nstfoutline!=-1 ){
        size_t pos = lines[nstxoutline].find("= ");
        size_t pos2 = lines[nstvoutline].find("= ");
        size_t pos3 = lines[nstfoutline].find("= ");
        size_t posnsteps = lines[nstepsline].find("= ");
        if (pos != std::string::npos && pos2 != std::string::npos && pos3 !=std::string::npos && nstepsline != -1 && posnsteps != std::string::npos){
            lines[nstxoutline] = lines[nstxoutline].substr(0,pos) + "= 1";
            lines[nstvoutline] = lines[nstvoutline].substr(0,pos2) + "= 1";
            lines[nstfoutline] = lines[nstfoutline].substr(0,pos3) + "= 1";
            lines[nstepsline] = lines[nstepsline].substr(0, posnsteps) + "= " + std::to_string(std::numeric_limits<int>::max());
            std::ofstream mdp("i.mdp");
            for (const auto&line:lines){
                mdp<<line<<std::endl;
            }
            mdp.close();
        }else{
            if (pos ==std::string::npos){
                std::cout<<"Fix gen-seed line (add = sign)"<<std::endl;
            }if (pos2 ==std::string::npos){
                std::cout<<"Fix dt line (add = sign)"<<std::endl;
            }if (pos3 ==std::string::npos){
                std::cout<<"Fix nsteps line (add = sign)" <<std::endl;
            }
        }
    }else{
        std::cout<<"insert nstxout, nstvout, nstfout lines into mdp file"<<std::endl;
    } 
}

void gromacsExport(Context* ctx,std::string groname, std::string cptfile, std::string outputname, int constanttopology){
    std::string topologyfile;
    std::string indexfile;
    if (constanttopology==1){
        topologyfile = "topol.top"; //constasnt topology requires this
        indexfile = "index.ndx";
    }else{
        topologyfile = groname.substr(0, groname.size()-4) + ".top";
        indexfile =  groname.substr(0, groname.size()-4) + ".ndx";
    }
    int rank =ctx->rank();
    if (rank ==0){
        InitializeGMX("i.mdp", groname, topologyfile, indexfile, "sil.tpr", cptfile, ctx);
        std::string cmd = "srun --cpu-bind=none -n 1 gmx_mpi mdrun -s sil.tpr -cpi " + cptfile + " -rerun " + cptfile + " -noappend -o sil.trr";
        std::system(cmd.c_str());
        std::vector<std::string> TRRS = listdirEnding(".", ".trr");
        std::string TRR = TRRS[0];
        FILE* process = popen(("srun  --overlap --overcommit --cpu-bind=none --mem-bind=none --mem=0 gmx_mpi trjconv -s sil.tpr -f " + TRR + " -o " + outputname + " -center -pbc mol").c_str(), "w");
        fputs("0\n0\n",process); //hopefully this works, in python it was a subprocess command that allowed multistep inputs
        pclose(process);
        DeleteExtraFiles(".", 3);
    }
    ctx->barrier();
}

void convertfromCPT(Context* ctx,std::string processedfolder, std::string rootdirectory, std::string outputfolder, int MDengine, int constanttopology){
    if (!std::filesystem::exists(outputfolder)){
        std::filesystem::create_directories(outputfolder);
    }
    //first create i.mdp
    if (MDengine ==1){
        generateimdp();
    }
    //okay the steps are to go through the cpt files in the newfolder, then for lammps we can instantly do the dump command. for gromacs find grofile etc.
    std::vector<std::string> cptfiles = listdirEnding(processedfolder, ".cpt");
    for (std::string file : cptfiles){
        std::string cptfile = processedfolder + "/" + file;
        int nom = std::stoi(file.substr(0,file.find('.')));
        if (MDengine ==0){
#ifdef USE_LAMMPS
            int argc = 0;
            char **argv = nullptr;
            LAMMPS_NS::LAMMPS* lmp = new LAMMPS_NS::LAMMPS(argc, argv, ctx->get_comm());
            lmp->input->one(("read_restart " + cptfile).c_str());
            lmp->input->one("reset_timestep 0");
            lmp->input->one( ("dump exportdump all custom 1 " + outputfolder + "/" + std::to_string(nom) + ".in id type x y z").c_str() );
            lmp->input->one("run 0");
            lmp->input->one("undump exportdump");
            delete lmp ;
#endif
            //lammps case, won't pass root directory in this case
        }else{
            std::string outputname = outputfolder + "/" + file.substr(0, file.find('.')) + ".gro";
            //gotta extract the int
            int gronumber = gronameparser(processedfolder, nom);
            std::string groname = rootdirectory + "/" + std::to_string(gronumber) + ".gro";
            gromacsExport(ctx,groname, cptfile, outputname, constanttopology);
        }
    }

}

int main(int argc, char** argv) { //okay here we get the inputs and decide the run type and set up the level of parallelism
    if (argc== 16){
        int MDenginetype = std::stoi(argv[1]);
        int simtype = std::stoi(argv[2]);
        std::string Outputfolder = argv[3];
        std::string runid = argv[4] ;
        double leftbound = std::stod(argv[5]);
        double rightbound = std::stod(argv[6]);
        std::string structurefiles = argv[7];
        std::string trail =argv[8];
        long long nsteps = std::stoll(argv[9]);
        std::string CUSTOMOP = argv[10];
        std::string rootdirectory = argv[11];
        int desiredtrajectories = std::stoi(argv[12]); //for probabilistic
        int strides = std::stoi(argv[13]);
        double deltaT = std::stod(argv[14]);
        int constanttopology = std::stoi(argv[15]);
        std::unique_ptr<Context> ctx; //initialize
        if (MDenginetype ==0){
#ifdef USE_LAMMPS
            MPI_Init(&argc, &argv);
            ctx = std::make_unique<ParallelContext>(MPI_COMM_WORLD);
#else
            std::cerr << "ERROR: Built without LAMMPS support" << std::endl;
            return 1;
#endif
        }else{
            //MPI_Init(&argc, &argv);
            ctx = std::make_unique<SerialContext>();
        }

        //std::this_thread::sleep_for(std::chrono::seconds(60)); //this is for requeue to avoid double running a process, I think this was causing requeue held, remove
        if (simtype ==0){
            Basin(ctx.get(), Outputfolder, runid, leftbound, rightbound, structurefiles, trail, nsteps, CUSTOMOP, rootdirectory, MDenginetype, strides,deltaT, constanttopology );
        }else if (simtype ==1){
            Probabilistic(ctx.get(), Outputfolder, runid, leftbound, rightbound, desiredtrajectories, structurefiles, trail, CUSTOMOP, rootdirectory, MDenginetype, strides,deltaT,constanttopology);
        }else if (simtype ==2){
            Exploration(ctx.get(), Outputfolder, runid, structurefiles, trail, nsteps, CUSTOMOP, rootdirectory, MDenginetype, strides,deltaT,constanttopology);
        }
    }else if( argc==1){
        Login();
    }
    else if(argc ==6){
        int simtype = std::stoi(argv[1]);
        std::string outputfolder = argv[2];
        double leftbound = std::stod(argv[3]);
        double rightbound = std::stod(argv[4]);
        int maxcrossings = std::stoi(argv[5]);
        controller(simtype, outputfolder, leftbound, rightbound, maxcrossings);
    }else if (argc==5){
        std::string RUNFOLDER = argv[1];
        std::string NEWFOLDER = argv[2];
        std::string prior = argv[3];
        int simtype = std::stoi(argv[4]);
        PostProcessing(RUNFOLDER, NEWFOLDER, prior, simtype);
    }else if (argc==7){
        std::unique_ptr<Context> ctx; //initialize
        MPI_Init(&argc, &argv);
        ctx = std::make_unique<SerialContext>(); // single runner
        std::string processedfolder = argv[2]; //skip first one it's meaningless to offset from controller
        std::string rootdirectory = argv[3];
        std::string outputfolder = argv[4];
        int MDEngine = std::stoi(argv[5]);
        int constanttopology = std::stoi(argv[6]);
        convertfromCPT(ctx.get(), processedfolder, rootdirectory, outputfolder, MDEngine, constanttopology);
    }
}
