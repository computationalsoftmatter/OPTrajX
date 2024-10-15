clear all; close all; clc

fid=fopen('0.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(1,i,:)=cell2mat(tMr); end
fid=fopen('1.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(2,i,:)=cell2mat(tMr); end
fid=fopen('2.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(3,i,:)=cell2mat(tMr); end
fid=fopen('3.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(4,i,:)=cell2mat(tMr); end
fid=fopen('4.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(5,i,:)=cell2mat(tMr); end
fid=fopen('5.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(6,i,:)=cell2mat(tMr); end
fid=fopen('6.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(7,i,:)=cell2mat(tMr); end
fid=fopen('7.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(8,i,:)=cell2mat(tMr); end
fid=fopen('8.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(9,i,:)=cell2mat(tMr); end

fid=fopen('9.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(10,i,:)=cell2mat(tMr); end
fid=fopen('10.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(11,i,:)=cell2mat(tMr); end
fid=fopen('11.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(12,i,:)=cell2mat(tMr); end
fid=fopen('12.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(13,i,:)=cell2mat(tMr); end
fid=fopen('13.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(14,i,:)=cell2mat(tMr); end
fid=fopen('14.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(15,i,:)=cell2mat(tMr); end
fid=fopen('15.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(16,i,:)=cell2mat(tMr); end
fid=fopen('16.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(17,i,:)=cell2mat(tMr); end
fid=fopen('17.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(18,i,:)=cell2mat(tMr); end
fid=fopen('18.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(19,i,:)=cell2mat(tMr); end

fid=fopen('19.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(20,i,:)=cell2mat(tMr); end
fid=fopen('20.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(21,i,:)=cell2mat(tMr); end
fid=fopen('21.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(22,i,:)=cell2mat(tMr); end
fid=fopen('22.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(23,i,:)=cell2mat(tMr); end
fid=fopen('23.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(24,i,:)=cell2mat(tMr); end
fid=fopen('24.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(25,i,:)=cell2mat(tMr); end
fid=fopen('25.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(26,i,:)=cell2mat(tMr); end
fid=fopen('26.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(27,i,:)=cell2mat(tMr); end
fid=fopen('27.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(28,i,:)=cell2mat(tMr); end
fid=fopen('28.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(29,i,:)=cell2mat(tMr); end

fid=fopen('29.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(30,i,:)=cell2mat(tMr); end
fid=fopen('30.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(31,i,:)=cell2mat(tMr); end
fid=fopen('31.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(32,i,:)=cell2mat(tMr); end
fid=fopen('32.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(33,i,:)=cell2mat(tMr); end
fid=fopen('33.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(34,i,:)=cell2mat(tMr); end
fid=fopen('34.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(35,i,:)=cell2mat(tMr); end
fid=fopen('35.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(36,i,:)=cell2mat(tMr); end
fid=fopen('36.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(37,i,:)=cell2mat(tMr); end
fid=fopen('37.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(38,i,:)=cell2mat(tMr); end
fid=fopen('38.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(39,i,:)=cell2mat(tMr); end

fid=fopen('39.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(40,i,:)=cell2mat(tMr); end
fid=fopen('40.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(41,i,:)=cell2mat(tMr); end
fid=fopen('41.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(42,i,:)=cell2mat(tMr); end
fid=fopen('42.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(43,i,:)=cell2mat(tMr); end
fid=fopen('43.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(44,i,:)=cell2mat(tMr); end
fid=fopen('44.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(45,i,:)=cell2mat(tMr); end
fid=fopen('45.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(46,i,:)=cell2mat(tMr); end
fid=fopen('46.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(47,i,:)=cell2mat(tMr); end
fid=fopen('47.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(48,i,:)=cell2mat(tMr); end
fid=fopen('48.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(49,i,:)=cell2mat(tMr); end

fid=fopen('49.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(50,i,:)=cell2mat(tMr); end
fid=fopen('50.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(51,i,:)=cell2mat(tMr); end
fid=fopen('51.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(52,i,:)=cell2mat(tMr); end
fid=fopen('52.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(53,i,:)=cell2mat(tMr); end
fid=fopen('53.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(54,i,:)=cell2mat(tMr); end
fid=fopen('54.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(55,i,:)=cell2mat(tMr); end
fid=fopen('55.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(56,i,:)=cell2mat(tMr); end
fid=fopen('56.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(57,i,:)=cell2mat(tMr); end
fid=fopen('57.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(58,i,:)=cell2mat(tMr); end
fid=fopen('58.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(59,i,:)=cell2mat(tMr); end

fid=fopen('59.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(60,i,:)=cell2mat(tMr); end
fid=fopen('60.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(61,i,:)=cell2mat(tMr); end
fid=fopen('61.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(62,i,:)=cell2mat(tMr); end
fid=fopen('62.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(63,i,:)=cell2mat(tMr); end
fid=fopen('63.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(64,i,:)=cell2mat(tMr); end
fid=fopen('64.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(65,i,:)=cell2mat(tMr); end
fid=fopen('65.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(66,i,:)=cell2mat(tMr); end
fid=fopen('66.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(67,i,:)=cell2mat(tMr); end
fid=fopen('67.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(68,i,:)=cell2mat(tMr); end
fid=fopen('68.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(69,i,:)=cell2mat(tMr); end

fid=fopen('69.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(70,i,:)=cell2mat(tMr); end
fid=fopen('70.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(71,i,:)=cell2mat(tMr); end
fid=fopen('71.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(72,i,:)=cell2mat(tMr); end
fid=fopen('72.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(73,i,:)=cell2mat(tMr); end
fid=fopen('73.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(74,i,:)=cell2mat(tMr); end
fid=fopen('74.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(75,i,:)=cell2mat(tMr); end
fid=fopen('75.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(76,i,:)=cell2mat(tMr); end
fid=fopen('76.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(77,i,:)=cell2mat(tMr); end
fid=fopen('77.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(78,i,:)=cell2mat(tMr); end
fid=fopen('78.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(79,i,:)=cell2mat(tMr); end

fid=fopen('79.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(80,i,:)=cell2mat(tMr); end
fid=fopen('80.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(81,i,:)=cell2mat(tMr); end
fid=fopen('81.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(82,i,:)=cell2mat(tMr); end
fid=fopen('82.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(83,i,:)=cell2mat(tMr); end
fid=fopen('83.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(84,i,:)=cell2mat(tMr); end
fid=fopen('84.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(85,i,:)=cell2mat(tMr); end
fid=fopen('85.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(86,i,:)=cell2mat(tMr); end
fid=fopen('86.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(87,i,:)=cell2mat(tMr); end
fid=fopen('87.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(88,i,:)=cell2mat(tMr); end
fid=fopen('88.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(89,i,:)=cell2mat(tMr); end

fid=fopen('89.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(90,i,:)=cell2mat(tMr); end
fid=fopen('90.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(91,i,:)=cell2mat(tMr); end
fid=fopen('91.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(92,i,:)=cell2mat(tMr); end
fid=fopen('92.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(93,i,:)=cell2mat(tMr); end
fid=fopen('93.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(94,i,:)=cell2mat(tMr); end
fid=fopen('94.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(95,i,:)=cell2mat(tMr); end
fid=fopen('95.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(96,i,:)=cell2mat(tMr); end
fid=fopen('96.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(97,i,:)=cell2mat(tMr); end
fid=fopen('97.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(98,i,:)=cell2mat(tMr); end
fid=fopen('98.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(99,i,:)=cell2mat(tMr); end

fid=fopen('99.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(100,i,:)=cell2mat(tMr);end
fid=fopen('100.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(101,i,:)=cell2mat(tMr); end
fid=fopen('101.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(102,i,:)=cell2mat(tMr); end
fid=fopen('102.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(103,i,:)=cell2mat(tMr); end
fid=fopen('103.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(104,i,:)=cell2mat(tMr); end
fid=fopen('104.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(105,i,:)=cell2mat(tMr); end
fid=fopen('105.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(106,i,:)=cell2mat(tMr); end
fid=fopen('106.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(107,i,:)=cell2mat(tMr); end
fid=fopen('107.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(108,i,:)=cell2mat(tMr); end
fid=fopen('108.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(109,i,:)=cell2mat(tMr); end

fid=fopen('109.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(110,i,:)=cell2mat(tMr); end
fid=fopen('110.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(111,i,:)=cell2mat(tMr); end
fid=fopen('111.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(112,i,:)=cell2mat(tMr); end
fid=fopen('112.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(113,i,:)=cell2mat(tMr); end
fid=fopen('113.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(114,i,:)=cell2mat(tMr); end
fid=fopen('114.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(115,i,:)=cell2mat(tMr); end
fid=fopen('115.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(116,i,:)=cell2mat(tMr); end
fid=fopen('116.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(117,i,:)=cell2mat(tMr); end
fid=fopen('117.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(118,i,:)=cell2mat(tMr); end
fid=fopen('118.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(119,i,:)=cell2mat(tMr); end

fid=fopen('119.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(120,i,:)=cell2mat(tMr); end
fid=fopen('120.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(121,i,:)=cell2mat(tMr); end
fid=fopen('121.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(122,i,:)=cell2mat(tMr); end
fid=fopen('122.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(123,i,:)=cell2mat(tMr); end
fid=fopen('123.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(124,i,:)=cell2mat(tMr); end
fid=fopen('124.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(125,i,:)=cell2mat(tMr); end
fid=fopen('125.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(126,i,:)=cell2mat(tMr); end
fid=fopen('126.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(127,i,:)=cell2mat(tMr); end
fid=fopen('127.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(128,i,:)=cell2mat(tMr); end
fid=fopen('128.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(129,i,:)=cell2mat(tMr); end

fid=fopen('129.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(130,i,:)=cell2mat(tMr); end
fid=fopen('130.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(131,i,:)=cell2mat(tMr); end
fid=fopen('131.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(132,i,:)=cell2mat(tMr); end
fid=fopen('132.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(133,i,:)=cell2mat(tMr); end
fid=fopen('133.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(134,i,:)=cell2mat(tMr); end
fid=fopen('134.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(135,i,:)=cell2mat(tMr); end
fid=fopen('135.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(136,i,:)=cell2mat(tMr); end
fid=fopen('136.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(137,i,:)=cell2mat(tMr); end
fid=fopen('137.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(138,i,:)=cell2mat(tMr); end
fid=fopen('138.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(139,i,:)=cell2mat(tMr); end

fid=fopen('139.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(140,i,:)=cell2mat(tMr); end
fid=fopen('140.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(141,i,:)=cell2mat(tMr); end
fid=fopen('141.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(142,i,:)=cell2mat(tMr); end
fid=fopen('142.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(143,i,:)=cell2mat(tMr); end
fid=fopen('143.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(144,i,:)=cell2mat(tMr); end
fid=fopen('144.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(145,i,:)=cell2mat(tMr); end
fid=fopen('145.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(146,i,:)=cell2mat(tMr); end
fid=fopen('146.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(147,i,:)=cell2mat(tMr); end
fid=fopen('147.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(148,i,:)=cell2mat(tMr); end
fid=fopen('148.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(149,i,:)=cell2mat(tMr); end

fid=fopen('149.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(150,i,:)=cell2mat(tMr); end
fid=fopen('150.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(151,i,:)=cell2mat(tMr); end
fid=fopen('151.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(152,i,:)=cell2mat(tMr); end
fid=fopen('152.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(153,i,:)=cell2mat(tMr); end
fid=fopen('153.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(154,i,:)=cell2mat(tMr); end
fid=fopen('154.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(155,i,:)=cell2mat(tMr); end
fid=fopen('155.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(156,i,:)=cell2mat(tMr); end
fid=fopen('156.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(157,i,:)=cell2mat(tMr); end
fid=fopen('157.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(158,i,:)=cell2mat(tMr); end
fid=fopen('158.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(159,i,:)=cell2mat(tMr); end

fid=fopen('159.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(160,i,:)=cell2mat(tMr); end
fid=fopen('160.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(161,i,:)=cell2mat(tMr); end
fid=fopen('161.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(162,i,:)=cell2mat(tMr); end
fid=fopen('162.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(163,i,:)=cell2mat(tMr); end
fid=fopen('163.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(164,i,:)=cell2mat(tMr); end
fid=fopen('164.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(165,i,:)=cell2mat(tMr); end
fid=fopen('165.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(166,i,:)=cell2mat(tMr); end
fid=fopen('166.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(167,i,:)=cell2mat(tMr); end
fid=fopen('167.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(168,i,:)=cell2mat(tMr); end
fid=fopen('168.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(169,i,:)=cell2mat(tMr); end

fid=fopen('169.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(170,i,:)=cell2mat(tMr); end
fid=fopen('170.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(171,i,:)=cell2mat(tMr); end
fid=fopen('171.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(172,i,:)=cell2mat(tMr); end
fid=fopen('172.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(173,i,:)=cell2mat(tMr); end
fid=fopen('173.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(174,i,:)=cell2mat(tMr); end
fid=fopen('174.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(175,i,:)=cell2mat(tMr); end
fid=fopen('175.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(176,i,:)=cell2mat(tMr); end
fid=fopen('176.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(177,i,:)=cell2mat(tMr); end
fid=fopen('177.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(178,i,:)=cell2mat(tMr); end
fid=fopen('178.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(179,i,:)=cell2mat(tMr); end

fid=fopen('179.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(180,i,:)=cell2mat(tMr); end
fid=fopen('180.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(181,i,:)=cell2mat(tMr); end
fid=fopen('181.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(182,i,:)=cell2mat(tMr); end
fid=fopen('182.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(183,i,:)=cell2mat(tMr); end
fid=fopen('183.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(184,i,:)=cell2mat(tMr); end
fid=fopen('184.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(185,i,:)=cell2mat(tMr); end
fid=fopen('185.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(186,i,:)=cell2mat(tMr); end
fid=fopen('186.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(187,i,:)=cell2mat(tMr); end
fid=fopen('187.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(188,i,:)=cell2mat(tMr); end
fid=fopen('188.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(189,i,:)=cell2mat(tMr); end

fid=fopen('189.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(190,i,:)=cell2mat(tMr); end
fid=fopen('190.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(191,i,:)=cell2mat(tMr); end
fid=fopen('191.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(192,i,:)=cell2mat(tMr); end
fid=fopen('192.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(193,i,:)=cell2mat(tMr); end
fid=fopen('193.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(194,i,:)=cell2mat(tMr); end
fid=fopen('194.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(195,i,:)=cell2mat(tMr); end
fid=fopen('195.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(196,i,:)=cell2mat(tMr); end
fid=fopen('196.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(197,i,:)=cell2mat(tMr); end
fid=fopen('197.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(198,i,:)=cell2mat(tMr); end
fid=fopen('198.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(199,i,:)=cell2mat(tMr); end

fid=fopen('199.pdb','r'); for i=1:497; tM=textscan(fid, '%s%f%s%s%f%f%f%f%f%f%f\r\n'); tMr=tM(6:8); M(200,i,:)=cell2mat(tMr); end

