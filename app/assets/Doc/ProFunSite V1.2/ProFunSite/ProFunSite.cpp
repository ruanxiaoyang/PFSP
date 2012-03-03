#pragma warning(disable:4018)
#pragma warning(disable:4996)
#include <process.h>
#include <windows.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <string>
#include <sstream>
#include <time.h>
#include <direct.h>
#include "sarray.h"
#include "darray.h"
using namespace std;
string TMPPATH;
bool BYSINPROGRESS=false;
#include "LogScore.h"
#include "CommonFunctions.h"
int CPUNUM=core_count();
#include "QuickSort.h"
#include "PAM.h"
#include "BLOSUM.h"
typedef sarray<char> chararray;
typedef darray<int> intmatrix;
typedef vector<intmatrix> vecintmatrix;
typedef sarray<int> intarray;
darray<int> Dayhoff250(getDayhoff250());
darray<double> Gonnet(GastonScorematx());
sarray<int> DayAccAAFreq(accumAAfreq('D'));
sarray<int> BloAccAAFreq(accumAAfreq('B'));
#include "display.h"
#include "ImpFile.h"
#include "EvoTree.h"
#include "GloAlign.h"
#include "RBLAST.h"
#include "BayesianDist.h"
#include "EvoDist.h"
#include "MultiAlign.h"
#include "exhaustwp.h"
#include "FunSite.h"

void main()
{
	HWND h=FindWindow("ConsoleWindowClass",NULL);
    SetWindowPos(h,HWND_TOP,0,0,1500,700,SWP_SHOWWINDOW);                           

	ifstream ChkLog("c:\\ProFunSite.log");                              //delete old log file
    if(ChkLog)
	{
		ChkLog.close();
	    system("del c:\\ProFunSite.log");
	}

	cout<<"         Protein Funtional Sites Prediction(PFSP) Version 1.2(Multi Core)"<<endl
		<<"              Thread safety powered by Intel(R) thread checker"<<endl
	    <<"                   Copyright belongs to Ruan Xiaoyang"<<endl
	    <<"                   Correspondence to ruansun@163.com"<<endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////Initiate scoring database
	cout<<"Initiating PAM Database....\n";
	vector<darray<double>> DayhoffPAMdb;
	darray<double> DayhoffPAM1(getDayhoffPAM1());
	CRPAMdb(DayhoffPAM1,DayhoffPAMdb,400);
	sarray<double> Dayhoffaafreq(getaafreq());
	vector<darray<double>> DayhoffROdb;                                 //relative odds database
	CRROdb(DayhoffPAMdb,Dayhoffaafreq,DayhoffROdb);
	vector<darray<double>> DayhoffScoredb;                              //Log(10) scoring database 
    CRScoredb(DayhoffROdb,10.0,DayhoffScoredb);
	cout<<"PAM Score Database Ready\n";

	cout<<"Initiating BLOSUM Database....\n";
	ifstream ChkBlkFile("c:\\windows\\blocks.dat");                     //block file
	if(!ChkBlkFile)
	{
		ifstream ImpBlkFile("blocks.dat");
		if(ImpBlkFile)
		   Copy("blocks.dat","c:\\windows\\blocks.dat");
		else
		{
			char blkpath[100];
impbf:		cout<<"Blocks file missing.Input path:";
			cin>>blkpath;
			ifstream ImpBlkFile(blkpath);
			if(!ImpBlkFile)
			{
				cout<<"Error!\n";
				goto impbf;
			}
			Copy(blkpath,"c:\\windows\\blocks.dat");
		}
		ImpBlkFile.close();
	}
	else
	    ChkBlkFile.close();

	ifstream ChkRFFile("c:\\windows\\_ProFunSit_BLOSUM_RltvFreq");
	if(!ChkRFFile)
	{
		ifstream ImpRFFile("_ProFunSit_BLOSUM_RltvFreq");
		if(ImpRFFile)
			Copy("_ProFunSit_BLOSUM_RltvFreq","c:\\windows\\_ProFunSit_BLOSUM_RltvFreq");
		else
		{
			char rfpath[100],skip;
    		cout<<"_ProFunSit_BLOSUM_RltvFreq file missing.Have you lost it?[Y/N]:";
			cin>>skip;
			if(skip=='N'||skip=='n')
			{
imprf:			cout<<"Input path:";
				cin>>rfpath;
			    ifstream ImpRFFile(rfpath);
			    if(!ImpRFFile)
			    {
					cout<<"Error!\n";
				    goto imprf;
				}
			    Copy(rfpath,"c:\\windows\\_ProFunSit_BLOSUM_RltvFreq");
			}
		}
		ImpRFFile.close();
	}
	else
	    ChkRFFile.close();
    cout<<"BLOSUM Score Database Ready\n";

    ifstream ChkTSFile("c:\\windows\\_ProFunSit_t_s_thld");
	if(!ChkTSFile)
	{
		ifstream ImpTSFile("_ProFunSit_t_s_thld");
		if(ImpTSFile)
			Copy("_ProFunSit_t_s_thld","c:\\windows\\_ProFunSit_t_s_thld");
		else
		{
			char TSpath[100],skip;
    		cout<<"_ProFunSit_t_s_thld file missing.Have you lost it?[Y/N]:";
			cin>>skip;
			if(skip=='N'||skip=='n')
			{
impts:			cout<<"Input path:";
				cin>>TSpath;
			    ifstream ImpTSFile(TSpath);
			    if(!ImpTSFile)
			    {cout<<"Error!\n";
				goto impts;}
			    Copy(TSpath,"c:\\windows\\_ProFunSit_t_s_thld");
			}
		}
		ImpTSFile.close();
	}
	else
	    ChkTSFile.close();
    readTSintomemory(string("c:\\windows\\_ProFunSit_t_s_thld"));

	ifstream ChkUSGD("c:\\windows\\PFSP_1.2_userguide.pdf");
	if(!ChkUSGD)
	{
		ifstream ImpUSGD("PFSP_1.2_userguide.pdf");
		if(ImpUSGD)
			system("copy PFSP_1.2_userguide.pdf c:\\windows\\PFSP_1.2_userguide.pdf");
        ImpUSGD.close();
	}
	else
	    ChkUSGD.close();

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////Initiating scoring database end
	cout<<"************************************************************\n"
		<<"NOTE:PFSP has two work modes available\n"
		<<"Fast mode:PFSP automatically handles most of the required parameters\n"
		<<"Custom mode:User will be asked to provide each parameter\n"
		<<"************************************************************\n";
mainmenu:   cout<<"Mainmenu\n"
				<<"Enter 1:Needleman global alignment\n"
				<<"      2:RBLAST glocal alignment\n"
				<<"      3:Construct custom BLOSUM scoring matrix from user specified block file\n"
				<<"      4:Multiple sequence alignment(MSA)\n"
				<<"      5:Protein funtional sites prediction\n"
				<<"      6:Protein funtional sites prediction(from MSA record)\n"
				<<"      I:Modify MSA result\n"
				<<"      H:A brief introduction to PFSP\n";
	cout<<"Input:";
	int gappnt,wlen,scmatxid,gamma=1,step;
	char choice, mode,format,fixscmatx,DorB,EXorAP,subtype;
	cin>>choice;
	proformat seqsandinf;
	char fpath[100];
	if(choice=='1')
	{
		cout<<"File format[F(fasta)/C(custom)]:";
		cin>>format;
ifpath0:cout<<"File path:";
        cin>>fpath;
		ifstream chkf(fpath);
		if(!chkf)
		{
			cout<<"No file found!\n";
			goto ifpath0;
		}
		if(format=='F'||format=='f')
           Fasta(fpath,seqsandinf);
		else if(format=='C'||format=='c')
           Custom(fpath,seqsandinf);
        cout<<"Fast/Custom mode[F/C]:";
        cin>>mode;
		if(mode=='F'||mode=='f')
			NWALN(seqsandinf,DayhoffScoredb,DayhoffScoredb[249],250,'N','P','D',-3);
        else if(mode=='C' ||mode=='c')
		{
			cout<<"Use fixed scoring matrix[Y/N]:";
			cin>>fixscmatx;
			if(fixscmatx=='Y'||fixscmatx=='y')
			{
				cout<<"Specify [scoring database][scoring matrix id][gap penalty] separate by space\n"
					<<"ep:D 250 -3 represent [PAM250][gap penalty -3]\n"
					<<"   B 62 -6 represent [BLOSUM62][gap penalty -6]\n"
					<<"Input:";
				cin>>DorB>>scmatxid>>gappnt;
				if(DorB=='D'||DorB=='d')
					NWALN(seqsandinf,DayhoffScoredb,DayhoffScoredb[scmatxid-1],scmatxid,'Y','S',DorB,gappnt);
				else
					NWALN(seqsandinf,DayhoffScoredb,BLOSUM(scmatxid,10),scmatxid,'Y','S',DorB,gappnt);
			}
			else
			{
				cout<<"Specify [approximate distance subtype][gap penalty] separate by space\n"
					<<"ep:P -3 represent [Poisson][gap penalty -3]\n"
					<<"   G -6 represent [Gamma][gap penalty -6]\n"
					<<"Input:";
				cin>>subtype>>gappnt;
				if(subtype=='G'||subtype=='g')
				{
					cout<<"Gamma distance engaged,specify gamma value:";
					cin>>gamma;
				}
				NWALN(seqsandinf,DayhoffScoredb,DayhoffScoredb[249],250,'N',subtype,'D',gappnt,gamma);
			}
		}
		chkf.close();
        seqsandinf._infvec.clear();seqsandinf._multiseqid.clear();
        goto mainmenu;
	}
	if(choice=='2')
	{
		cout<<"File format[F(fasta)/C(custom)]:";
        cin>>format;
ifpath1:cout<<"File path:";
        cin>>fpath;
		ifstream chkf(fpath);
		if(!chkf)
		{
			cout<<"No file found!\n";
			goto ifpath1;
		}
		if(format=='F'||format=='f')
           Fasta(fpath,seqsandinf);
		else if(format=='C'||format=='c')
           Custom(fpath,seqsandinf);
        cout<<"Fast/Custom mode[F/C]:";
        cin>>mode;
		if(mode=='F'||mode=='f')
			RBLAST(seqsandinf,DayhoffScoredb,DayhoffScoredb[249],250,'N','D','A','P',-3);
		else if(mode=='C' ||mode=='c')
		{
			cout<<"Use fixed scoring matrix[Y/N]:";
			cin>>fixscmatx;
			if(fixscmatx=='Y'||fixscmatx=='y')
			{
				cout<<"Specify [scoring database][matrix id][gap penalty][word length] separate by space\n"
					<<"ep:D 250 -3 3 represent [PAM250][gap penalty -3][word length 3]\n"
					<<"  :B 62 -6 3 represent [BLOSUM62][gap penalty -6][word length 3]\n"
					<<"Input:";
				cin>>DorB>>scmatxid>>gappnt>>wlen;
				if(DorB=='D'||DorB=='d')
					RBLAST(seqsandinf,DayhoffScoredb,DayhoffScoredb[scmatxid-1],scmatxid,'Y','D','E','B',gappnt,1,wlen);
				else
					RBLAST(seqsandinf,DayhoffScoredb,BLOSUM(scmatxid,10),scmatxid,'Y','B','E','B',gappnt,1,wlen);
			}
			else if(fixscmatx=='N'||fixscmatx=='n')
			{
				cout<<"Specify [distance estimation type][subtype][gap penalty][word length] separate by space\n"
					<<"ep:E B -3 3 represent [Bayesian exact][BLAST local][gap penalty -3][word length 3]\n"
					<<"ep:A P -6 3 represent [approximate][Poisson][gap penalty -6][word length 3]\n"
					<<"Input:";
				cin>>EXorAP>>subtype>>gappnt>>wlen;
                if(subtype=='G'||subtype=='g')
				{
					cout<<"Gamma distance engaged,specify gamma value:";
					cin>>gamma;
				}
				RBLAST(seqsandinf,DayhoffScoredb,DayhoffScoredb[249],250,'N','D',EXorAP,subtype,gappnt,gamma,wlen);
			}
		}
		chkf.close();
		seqsandinf._infvec.clear();seqsandinf._multiseqid.clear();
		goto mainmenu;
	}
	int cvgstd,logbase;
	if(choice=='3')
	{
		cout<<"Specify [reclustering standard(100~30)][log base] separate by space\n"
			<<"ep:62 10 represent [reclustering standard 62][log base 10]\n"
			<<"Input:";
		cin>>cvgstd>>logbase;
		BLOSUM(cvgstd,logbase,'N');
		goto mainmenu;
	}
	int seqnumthld;
	darray<char> dspmatx;
	darray<int> seqsidmatx;
	sarray<char> identical;
	if(choice=='4' ||choice=='5')
	{
        cout<<"File format[F(fasta)/C(custom)]:";
        cin>>format;
ifpath2:cout<<"File path:";
        cin>>fpath;
		ifstream chkf(fpath);
		if(!chkf)
		{
			cout<<"No file found!\n";
			goto ifpath2;
		}
		if(format=='F'||format=='f')
           Fasta(fpath,seqsandinf);
		else if(format=='C'||format=='c')
           Custom(fpath,seqsandinf);
		if(seqsandinf._seqnum<3)
		{
			cout<<"At least 3 sequences are needed to enable MSA\n";
			goto mainmenu;
		}
		if(choice=='5' && seqsandinf._seqnum<20)
		{
			cout<<"At least 20 sequences are needed to enable funtional site prediction\n";
			goto mainmenu;
		}
		cout<<"Fast/Custom mode[F/C]:";
        cin>>mode;
		if(mode=='F'||mode=='f')
		   MSA(seqsandinf,DayhoffScoredb,seqsidmatx,dspmatx,identical,'A','P',3,-3,50,2);
		else if(mode=='C'||mode=='c')
		{
		   cout<<"Specify [distance estimation type][subtype][gap penalty][word length] separate by space\n"
			   <<"ep:E B -3 3 represent [Bayesian exact][BLAST local][gap penalty -3][word length 3]\n"
			   <<"ep:A P -6 3 represent [Approximate][Poisson][gap penalty -6][word length 3]\n"
			   <<"Input:";
		   cin>>EXorAP>>subtype>>gappnt>>wlen;
		   if(EXorAP=='E'||EXorAP=='e')
		   {
			   cout<<"Specify bayesian exact distance step(large step means higher speed but lower accuracy,recommend 10):";
			   cin>>step;
		   }
		   if(subtype=='G'||subtype=='g')
			{
				cout<<"Gamma distance engaged,specify gamma value:";
				cin>>gamma;
			}
		    MSA(seqsandinf,DayhoffScoredb,seqsidmatx,dspmatx,identical,EXorAP,subtype,wlen,gappnt,50,2,95,'D',gamma,step,4);
		}
		if(choice=='5')
		{
			cout<<"Specify sequence number threshold(>=10):";
			cin>>seqnumthld;
			PFSP(seqsidmatx,dspmatx,identical,seqnumthld);
		}
		dspmatx.clear();
		seqsidmatx.clear();
		identical.clear();
		chkf.close();
		goto mainmenu;
	}	   
	if(choice=='6')
	{
ifpath4:cout<<"File path of MSA result:";
        cin>>fpath;
		ifstream chkf(fpath);
		if(!chkf)
		{
			cout<<"No file found!\n";
			goto ifpath4;
		}
		cout<<"Specify sequence number threshold(>=10):";
		cin>>seqnumthld;
		if(PFSP(fpath,seqnumthld)==0)
			cout<<"At leat 20 sequences are needed to enable funtional site prediction\n";
		chkf.close();
		goto mainmenu;
	}
    if(choice=='I' || choice=='i')
	{
ifpath5:cout<<"File path of MSA result:";
        cin>>fpath;
		ifstream chkf(fpath);
		if(!chkf)
		{
			cout<<"No file found!\n";
			goto ifpath5;
		}
		_ImpMSA(fpath,dspmatx,identical);
        _MdfMSA(dspmatx);
		goto mainmenu;
	}
	if(choice=='H'|| choice=='h')
	{
		cout<<"See attaced PDF instruction file for detailed information.Open PDF instruction file?[Y/N]:";
		char answer;
		cin>>answer;
		if(answer=='Y'||answer=='y')
		{
			TMPPATH="c:\\windows\\PFSP_1.2_userguide.pdf";
			DSPFILE(TMPPATH);
		}
		goto mainmenu;
	}
}