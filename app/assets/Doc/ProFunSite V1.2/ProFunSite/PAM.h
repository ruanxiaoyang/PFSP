//////////////////////////////////////////////////////////////
//           Point Accepted Mutation (PAM)                  //
//        M.O. Dayhoff, R.M. Schwartz, and B. C, Orcutt     //
//      1978 A model of evolutionary change in proteins     //
//           Copyright Belongs to Ruan Xiaoyang             //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
//This module intended to calculate PAM score matrix from PAM1 to PAM400.A log transformation function is available to transform relative odds into logN based score
//The most useful output is contained in the CRScoredb function,which output vector<darray<double>> score database in the 3rd parameter  
//A log file was saved at c:\ProFunSite.log

darray<double> getDayhoffPAM1()           //PAM1 frequency data from Dayhoff literature
{
      darray<double> PAM1(20,20,0.0);
//A                R                N                D                C                Q                E                G                H                I                L                 K                 M                 F                 P                 S                 T                 W                 Y                 V 
PAM1(0,0)=0.9867;PAM1(0,1)=0.0002;PAM1(0,2)=0.0009;PAM1(0,3)=0.0010;PAM1(0,4)=0.0003;PAM1(0,5)=0.0008;PAM1(0,6)=0.0017;PAM1(0,7)=0.0021;PAM1(0,8)=0.0002;PAM1(0,9)=0.0006;PAM1(0,10)=0.0004;PAM1(0,11)=0.0002;PAM1(0,12)=0.0006;PAM1(0,13)=0.0002;PAM1(0,14)=0.0022;PAM1(0,15)=0.0035;PAM1(0,16)=0.0032;PAM1(0,17)=0.0000;PAM1(0,18)=0.0002;PAM1(0,19)=0.0018;
PAM1(1,0)=0.0001;PAM1(1,1)=0.9913;PAM1(1,2)=0.0001;PAM1(1,3)=0.0000;PAM1(1,4)=0.0001;PAM1(1,5)=0.0010;PAM1(1,6)=0.0000;PAM1(1,7)=0.0000;PAM1(1,8)=0.0010;PAM1(1,9)=0.0003;PAM1(1,10)=0.0001;PAM1(1,11)=0.0019;PAM1(1,12)=0.0004;PAM1(1,13)=0.0001;PAM1(1,14)=0.0004;PAM1(1,15)=0.0006;PAM1(1,16)=0.0001;PAM1(1,17)=0.0008;PAM1(1,18)=0.0000;PAM1(1,19)=0.0001;
PAM1(2,0)=0.0004;PAM1(2,1)=0.0001;PAM1(2,2)=0.9822;PAM1(2,3)=0.0036;PAM1(2,4)=0.0000;PAM1(2,5)=0.0004;PAM1(2,6)=0.0006;PAM1(2,7)=0.0006;PAM1(2,8)=0.0021;PAM1(2,9)=0.0003;PAM1(2,10)=0.0001;PAM1(2,11)=0.0013;PAM1(2,12)=0.0000;PAM1(2,13)=0.0001;PAM1(2,14)=0.0002;PAM1(2,15)=0.0020;PAM1(2,16)=0.0009;PAM1(2,17)=0.0001;PAM1(2,18)=0.0004;PAM1(2,19)=0.0001;
PAM1(3,0)=0.0006;PAM1(3,1)=0.0000;PAM1(3,2)=0.0042;PAM1(3,3)=0.9859;PAM1(3,4)=0.0000;PAM1(3,5)=0.0006;PAM1(3,6)=0.0053;PAM1(3,7)=0.0006;PAM1(3,8)=0.0004;PAM1(3,9)=0.0001;PAM1(3,10)=0.0000;PAM1(3,11)=0.0003;PAM1(3,12)=0.0000;PAM1(3,13)=0.0000;PAM1(3,14)=0.0001;PAM1(3,15)=0.0005;PAM1(3,16)=0.0003;PAM1(3,17)=0.0000;PAM1(3,18)=0.0000;PAM1(3,19)=0.0001;
PAM1(4,0)=0.0001;PAM1(4,1)=0.0001;PAM1(4,2)=0.0000;PAM1(4,3)=0.0000;PAM1(4,4)=0.9973;PAM1(4,5)=0.0000;PAM1(4,6)=0.0000;PAM1(4,7)=0.0000;PAM1(4,8)=0.0001;PAM1(4,9)=0.0001;PAM1(4,10)=0.0000;PAM1(4,11)=0.0000;PAM1(4,12)=0.0000;PAM1(4,13)=0.0000;PAM1(4,14)=0.0001;PAM1(4,15)=0.0005;PAM1(4,16)=0.0001;PAM1(4,17)=0.0000;PAM1(4,18)=0.0003;PAM1(4,19)=0.0002;
PAM1(5,0)=0.0003;PAM1(5,1)=0.0009;PAM1(5,2)=0.0004;PAM1(5,3)=0.0005;PAM1(5,4)=0.0000;PAM1(5,5)=0.9876;PAM1(5,6)=0.0027;PAM1(5,7)=0.0001;PAM1(5,8)=0.0023;PAM1(5,9)=0.0001;PAM1(5,10)=0.0003;PAM1(5,11)=0.0006;PAM1(5,12)=0.0004;PAM1(5,13)=0.0000;PAM1(5,14)=0.0006;PAM1(5,15)=0.0002;PAM1(5,16)=0.0002;PAM1(5,17)=0.0000;PAM1(5,18)=0.0000;PAM1(5,19)=0.0001;
PAM1(6,0)=0.0010;PAM1(6,1)=0.0000;PAM1(6,2)=0.0007;PAM1(6,3)=0.0056;PAM1(6,4)=0.0000;PAM1(6,5)=0.0035;PAM1(6,6)=0.9865;PAM1(6,7)=0.0004;PAM1(6,8)=0.0002;PAM1(6,9)=0.0003;PAM1(6,10)=0.0001;PAM1(6,11)=0.0004;PAM1(6,12)=0.0001;PAM1(6,13)=0.0000;PAM1(6,14)=0.0003;PAM1(6,15)=0.0004;PAM1(6,16)=0.0002;PAM1(6,17)=0.0000;PAM1(6,18)=0.0001;PAM1(6,19)=0.0002;
PAM1(7,0)=0.0021;PAM1(7,1)=0.0001;PAM1(7,2)=0.0012;PAM1(7,3)=0.0011;PAM1(7,4)=0.0001;PAM1(7,5)=0.0003;PAM1(7,6)=0.0007;PAM1(7,7)=0.9935;PAM1(7,8)=0.0001;PAM1(7,9)=0.0000;PAM1(7,10)=0.0001;PAM1(7,11)=0.0002;PAM1(7,12)=0.0001;PAM1(7,13)=0.0001;PAM1(7,14)=0.0003;PAM1(7,15)=0.0021;PAM1(7,16)=0.0003;PAM1(7,17)=0.0000;PAM1(7,18)=0.0000;PAM1(7,19)=0.0005;
PAM1(8,0)=0.0001;PAM1(8,1)=0.0008;PAM1(8,2)=0.0018;PAM1(8,3)=0.0003;PAM1(8,4)=0.0001;PAM1(8,5)=0.0020;PAM1(8,6)=0.0001;PAM1(8,7)=0.0000;PAM1(8,8)=0.9912;PAM1(8,9)=0.0000;PAM1(8,10)=0.0001;PAM1(8,11)=0.0001;PAM1(8,12)=0.0000;PAM1(8,13)=0.0002;PAM1(8,14)=0.0003;PAM1(8,15)=0.0001;PAM1(8,16)=0.0001;PAM1(8,17)=0.0001;PAM1(8,18)=0.0004;PAM1(8,19)=0.0001;
PAM1(9,0)=0.0002;PAM1(9,1)=0.0002;PAM1(9,2)=0.0003;PAM1(9,3)=0.0001;PAM1(9,4)=0.0002;PAM1(9,5)=0.0001;PAM1(9,6)=0.0002;PAM1(9,7)=0.0000;PAM1(9,8)=0.0000;PAM1(9,9)=0.9872;PAM1(9,10)=0.0009;PAM1(9,11)=0.0002;PAM1(9,12)=0.0012;PAM1(9,13)=0.0007;PAM1(9,14)=0.0000;PAM1(9,15)=0.0001;PAM1(9,16)=0.0007;PAM1(9,17)=0.0000;PAM1(9,18)=0.0001;PAM1(9,19)=0.0033;
PAM1(10,0)=0.0003;PAM1(10,1)=0.0001;PAM1(10,2)=0.0003;PAM1(10,3)=0.0000;PAM1(10,4)=0.0000;PAM1(10,5)=0.0006;PAM1(10,6)=0.0001;PAM1(10,7)=0.0001;PAM1(10,8)=0.0004;PAM1(10,9)=0.0022;PAM1(10,10)=0.9947;PAM1(10,11)=0.0002;PAM1(10,12)=0.0045;PAM1(10,13)=0.0013;PAM1(10,14)=0.0003;PAM1(10,15)=0.0001;PAM1(10,16)=0.0003;PAM1(10,17)=0.0004;PAM1(10,18)=0.0002;PAM1(10,19)=0.0015;
PAM1(11,0)=0.0002;PAM1(11,1)=0.0037;PAM1(11,2)=0.0025;PAM1(11,3)=0.0006;PAM1(11,4)=0.0000;PAM1(11,5)=0.0012;PAM1(11,6)=0.0007;PAM1(11,7)=0.0002;PAM1(11,8)=0.0002;PAM1(11,9)=0.0004;PAM1(11,10)=0.0001;PAM1(11,11)=0.9926;PAM1(11,12)=0.0020;PAM1(11,13)=0.0000;PAM1(11,14)=0.0003;PAM1(11,15)=0.0008;PAM1(11,16)=0.0011;PAM1(11,17)=0.0000;PAM1(11,18)=0.0001;PAM1(11,19)=0.0001;
PAM1(12,0)=0.0001;PAM1(12,1)=0.0001;PAM1(12,2)=0.0000;PAM1(12,3)=0.0000;PAM1(12,4)=0.0000;PAM1(12,5)=0.0002;PAM1(12,6)=0.0000;PAM1(12,7)=0.0000;PAM1(12,8)=0.0000;PAM1(12,9)=0.0005;PAM1(12,10)=0.0008;PAM1(12,11)=0.0004;PAM1(12,12)=0.9874;PAM1(12,13)=0.0001;PAM1(12,14)=0.0000;PAM1(12,15)=0.0001;PAM1(12,16)=0.0002;PAM1(12,17)=0.0000;PAM1(12,18)=0.0000;PAM1(12,19)=0.0004;
PAM1(13,0)=0.0001;PAM1(13,1)=0.0001;PAM1(13,2)=0.0001;PAM1(13,3)=0.0000;PAM1(13,4)=0.0000;PAM1(13,5)=0.0000;PAM1(13,6)=0.0000;PAM1(13,7)=0.0001;PAM1(13,8)=0.0002;PAM1(13,9)=0.0008;PAM1(13,10)=0.0006;PAM1(13,11)=0.0000;PAM1(13,12)=0.0004;PAM1(13,13)=0.9946;PAM1(13,14)=0.0000;PAM1(13,15)=0.0002;PAM1(13,16)=0.0001;PAM1(13,17)=0.0003;PAM1(13,18)=0.0028;PAM1(13,19)=0.0000;
PAM1(14,0)=0.0013;PAM1(14,1)=0.0005;PAM1(14,2)=0.0002;PAM1(14,3)=0.0001;PAM1(14,4)=0.0001;PAM1(14,5)=0.0008;PAM1(14,6)=0.0003;PAM1(14,7)=0.0002;PAM1(14,8)=0.0005;PAM1(14,9)=0.0001;PAM1(14,10)=0.0002;PAM1(14,11)=0.0002;PAM1(14,12)=0.0001;PAM1(14,13)=0.0001;PAM1(14,14)=0.9926;PAM1(14,15)=0.0012;PAM1(14,16)=0.0004;PAM1(14,17)=0.0000;PAM1(14,18)=0.0000;PAM1(14,19)=0.0002;
PAM1(15,0)=0.0028;PAM1(15,1)=0.0011;PAM1(15,2)=0.0034;PAM1(15,3)=0.0007;PAM1(15,4)=0.0011;PAM1(15,5)=0.0004;PAM1(15,6)=0.0006;PAM1(15,7)=0.0016;PAM1(15,8)=0.0002;PAM1(15,9)=0.0002;PAM1(15,10)=0.0001;PAM1(15,11)=0.0007;PAM1(15,12)=0.0004;PAM1(15,13)=0.0003;PAM1(15,14)=0.0017;PAM1(15,15)=0.9840;PAM1(15,16)=0.0038;PAM1(15,17)=0.0005;PAM1(15,18)=0.0002;PAM1(15,19)=0.0002;
PAM1(16,0)=0.0022;PAM1(16,1)=0.0002;PAM1(16,2)=0.0013;PAM1(16,3)=0.0004;PAM1(16,4)=0.0001;PAM1(16,5)=0.0003;PAM1(16,6)=0.0002;PAM1(16,7)=0.0002;PAM1(16,8)=0.0001;PAM1(16,9)=0.0011;PAM1(16,10)=0.0002;PAM1(16,11)=0.0008;PAM1(16,12)=0.0006;PAM1(16,13)=0.0001;PAM1(16,14)=0.0005;PAM1(16,15)=0.0032;PAM1(16,16)=0.9871;PAM1(16,17)=0.0000;PAM1(16,18)=0.0002;PAM1(16,19)=0.0009;
PAM1(17,0)=0.0000;PAM1(17,1)=0.0002;PAM1(17,2)=0.0000;PAM1(17,3)=0.0000;PAM1(17,4)=0.0000;PAM1(17,5)=0.0000;PAM1(17,6)=0.0000;PAM1(17,7)=0.0000;PAM1(17,8)=0.0000;PAM1(17,9)=0.0000;PAM1(17,10)=0.0000;PAM1(17,11)=0.0000;PAM1(17,12)=0.0000;PAM1(17,13)=0.0001;PAM1(17,14)=0.0000;PAM1(17,15)=0.0001;PAM1(17,16)=0.0000;PAM1(17,17)=0.9976;PAM1(17,18)=0.0001;PAM1(17,19)=0.0000;
PAM1(18,0)=0.0001;PAM1(18,1)=0.0000;PAM1(18,2)=0.0003;PAM1(18,3)=0.0000;PAM1(18,4)=0.0003;PAM1(18,5)=0.0000;PAM1(18,6)=0.0001;PAM1(18,7)=0.0000;PAM1(18,8)=0.0004;PAM1(18,9)=0.0001;PAM1(18,10)=0.0001;PAM1(18,11)=0.0000;PAM1(18,12)=0.0000;PAM1(18,13)=0.0021;PAM1(18,14)=0.0000;PAM1(18,15)=0.0001;PAM1(18,16)=0.0001;PAM1(18,17)=0.0002;PAM1(18,18)=0.9945;PAM1(18,19)=0.0001;
PAM1(19,0)=0.0013;PAM1(19,1)=0.0002;PAM1(19,2)=0.0001;PAM1(19,3)=0.0001;PAM1(19,4)=0.0003;PAM1(19,5)=0.0002;PAM1(19,6)=0.0002;PAM1(19,7)=0.0003;PAM1(19,8)=0.0003;PAM1(19,9)=0.0057;PAM1(19,10)=0.0011;PAM1(19,11)=0.0001;PAM1(19,12)=0.0017;PAM1(19,13)=0.0001;PAM1(19,14)=0.0003;PAM1(19,15)=0.0002;PAM1(19,16)=0.0010;PAM1(19,17)=0.0000;PAM1(19,18)=0.0002;PAM1(19,19)=0.9901;
return PAM1;
}


sarray<double> getaafreq(const char & _choice='D')                   //AA frequency data from Dayhoff/Blosum100
{
    sarray<double> aafreq(20,0.0);
    if(_choice=='D'|| _choice=='d')
	{
//A                R                N                D                C                Q                E                G                H                I                L                 K                 M                 F                 P                 S                 T                 W                 Y                 V 
aafreq[0]=0.087; aafreq[1]=0.041; aafreq[2]=0.04;  aafreq[3]=0.047; aafreq[4]=0.033; aafreq[5]=0.038; aafreq[6]=0.05;  aafreq[7]=0.089; aafreq[8]=0.034; aafreq[9]=0.037; aafreq[10]=0.085; aafreq[11]=0.081; aafreq[12]=0.015; aafreq[13]=0.04;  aafreq[14]=0.051; aafreq[15]=0.07;  aafreq[16]=0.058; aafreq[17]=0.01;  aafreq[18]=0.03;  aafreq[19]=0.065;
	}
	else if(_choice=='B' || _choice=='b')
	{
aafreq[0]=0.073; aafreq[1]=0.049; aafreq[2]=0.042; aafreq[3]=0.054; aafreq[4]=0.029; aafreq[5]=0.034; aafreq[6]=0.057; aafreq[7]=0.081; aafreq[8]=0.025; aafreq[9]=0.065; aafreq[10]=0.095; aafreq[11]=0.052; aafreq[12]=0.025; aafreq[13]=0.046; aafreq[14]=0.039; aafreq[15]=0.058; aafreq[16]=0.053; aafreq[17]=0.015; aafreq[18]=0.035; aafreq[19]=0.073;
	}
    return aafreq;
}
sarray<int> accumAAfreq(char _choice='D')    //Dayhoff/Blosum100 AA accumulate frequency data.Used for creating random AAid sequence
{
	sarray<int> accumfreq(20,0);
	if(_choice=='D' || _choice=='d'||_choice=='G'||_choice=='g')
	{
	    accumfreq[0]=87;accumfreq[1]=128;accumfreq[2]=168;accumfreq[3]=215;accumfreq[4]=248;accumfreq[5]=286;accumfreq[6]=336;accumfreq[7]=425;accumfreq[8]=459;accumfreq[9]=496;accumfreq[10]=581;
	    accumfreq[11]=662;accumfreq[12]=677;accumfreq[13]=717;accumfreq[14]=768;accumfreq[15]=838;accumfreq[16]=896;accumfreq[17]=906;accumfreq[18]=936;accumfreq[19]=1000;
	}
	else if(_choice='B' || _choice=='b')
	{
		accumfreq[0]=73;accumfreq[1]=122;accumfreq[2]=164;accumfreq[3]=218;accumfreq[4]=247;accumfreq[5]=281;accumfreq[6]=338;accumfreq[7]=419;accumfreq[8]=444;accumfreq[9]=509;accumfreq[10]=604;
		accumfreq[11]=656;accumfreq[12]=681;accumfreq[13]=727;accumfreq[14]=766;accumfreq[15]=824;accumfreq[16]=877;accumfreq[17]=892;accumfreq[18]=927;accumfreq[19]=1000;
	}
        return accumfreq;
}
//generate random AA id sequence.optimized for 20 kinds of AA.it costs averagely (n^2+n+20)/2n time to generate one AA.when n=4 it has the smallest value (5 computer clock) 
void randAAidseq(const int & _length,sarray<int> & _accumdb,sarray<int> & randseqid)    
{
	int temp;
	randseqid.resize(_length);
    for(int n=0;n<_length;++n)
	{
		temp=rand()%1001;
		if(temp<=_accumdb[4])
		{
			for(int i=0;i<=4;++i)
	        {
		       if(temp<=_accumdb[i])
		       {
				   randseqid[n]=i;
			      break;
		       }
	        }
		}
		else if(temp<=_accumdb[9])
		{
			for(int i=5;i<=9;++i)
			{
			   if(temp<=_accumdb[i])
		       {
			      randseqid[n]=i;
			      break;
		       }
			}
		}
		else if(temp<=_accumdb[14])
		{
			for(int i=10;i<=14;++i)
			{
			   if(temp<=_accumdb[i])
		       {
			      randseqid[n]=i;
			      break;
		       }
			}
		}
		else
		{
			for(int i=15;i<=19;++i)
			{
			   if(temp<=_accumdb[i])
		       {
			      randseqid[n]=i;
			      break;
		       }
			}
		}
	}
}

void CRPAMdb(darray<double> & _PAM1,vector<darray<double> > & _PAMdb,const int & _updist)      //create mutation frequency database from PAM 1 to 400
{
	_PAMdb.push_back(_PAM1);
	for(int i=1;i<_updist;++i)
		_PAMdb.push_back(_PAMdb[i-1].multiply(_PAMdb[0]));                 //matrix multiplication
}

void CRROdb(vector<darray<double>> & _PAMdb,sarray<double> & _aafreq,vector<darray<double>> & _RelatedOddsdb)    //create related odds database
{
    darray<double> ROtemp(20,20);
    for(int s=0;s<_PAMdb.size();++s)
	{
		for(int i=0;i<20;++i)
		{
			for(int j=i;j<20;++j)
			{
                if(i==j)
	            ROtemp(i,j)=_PAMdb[s](i,j)/_aafreq[i];
				else
				{
					ROtemp(i,j)=pow(_PAMdb[s](i,j)*_PAMdb[s](j,i)/(_aafreq[i]*_aafreq[j]),0.5);
					ROtemp(j,i)=ROtemp(i,j);
				}
			}
		}
		_RelatedOddsdb.push_back(ROtemp);
		ROtemp.fast_resize(20,20);
	}
}
void CRScoredb(vector<darray<double>> & _RelatedOddsdb,const double & _base,vector<darray<double>> & _Scoredb)      //create score database
{
	clock_t start,end;
	start=clock();
	darray<double> scoretemp(20,20,0);
	for(int s=0;s<_RelatedOddsdb.size();++s)
	{
		for(int i=0;i<20;++i)
		{
			for(int j=i;j<20;++j)
			{
                if(_RelatedOddsdb[s](i,j)>0)
				   scoretemp(i,j)=10*log(_RelatedOddsdb[s](i,j))/log(_base);
				else
                   scoretemp(i,j)=-50;
				if(i!=j)
				   scoretemp(j,i)=scoretemp(i,j);
			}
		}
		_Scoredb.push_back(scoretemp);
		scoretemp.fast_resize(20,20);
	}
	end=clock();
	char time[10],date[10];
	_strtime(time);_strdate(date);
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<" Process time "<<end-start<<" minisec"<<endl;
	outlog<<"log "<<_base<<" PAM score database (1-"<<_RelatedOddsdb.size()<<") created"<<endl<<endl;
}


