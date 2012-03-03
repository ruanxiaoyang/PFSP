//////////////////////////////////////////////////////////////
//             Coupling Energy Estimation                   //
//              Steve W. Lockless, et al.                   //
//               Science 286, 295 (1999);                   //
//           Copyright Belongs to Ruan Xiaoyang             //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
//This module estimate coupling energy for protein sequence.It calculates both detag(estimated conservative degree) and deta_detag(estimate detag change if other site mutated) 
//For each run,a log file was saved at c:\ProFunSite.log indicating the sequence number of the processed MSA

//Menu:[PFSP]-[_LogRLPmatx]-[LogBinoP]
//           -[_CalDetaG]
//           -[QuickSort]
//           -[_ImpMSA]

struct funsitedata                                                              //struct hold all deta_detaG information of a MSA
{
	sarray<int> _funsitearr;                                                    //store id of site which have mutation
	vector<sarray<char>> _AAvec;                                                //sarray<char> contain AA mutations at each site
	vector<sarray<double>> _avgD_DGvec;
	vector<darray<double>> _D_DGvec;                                            //darray<double> contain coupling energy of the rest sites when one site mutated
};
//Menu:[PFSP]-[_ImpMSA]
//Import MSA file from existing data and store into "dspmatx".
//Take 1 value."_path"the location of MSA result
//Output 2 values."dspmatx"display matrix,aligned AAid seqs."identical"site with identical AA
void _ImpMSA(char * _path,darray<char> & dspmatx,sarray<char> & identical)
{
	char a;
	int row,column,i=0,j=0;
	ifstream imp(_path);
	imp>>row>>column;                                                           //the first row of MSA file store row and column number
	identical.resize(column);
	dspmatx.fast_resize(row,column);
	imp.get(a);
	while(imp.get(a))
	{
		if(a=='\n')
			break;
		else
		{
			identical[j]=a;
			j+=1;
		}
	}
	j=0;
	while(imp.get(a))
	{
		if(a=='\n')
		{
			i+=1;
			j=0;
		}
		else
		{
			dspmatx(i,j)=a;
			j+=1;
		}
	}
	imp.close();
}
//Menu:[PFSP]-[_ImpMSAClustalW]
//Import from ClustalW generated MSA file and store into "dspmatx".
//Take 1 value."_path"the location of MSA result
//Output 2 values."dspmatx"display matrix,aligned AAid seqs."identical"site with identical AA
void _ImpMSAClustalW(char * _path,darray<char> & dspmatx,sarray<char> & identical)
{
     char a;
	 string x;
	 bool enter=false;
	 int rownum=0;
	 ifstream imp(_path);
	 getline(imp,x);
	 while(imp.get(a))
	 {
		 if(a!=' ' && a!='\n')
		 {
			 imp.unget();
     		 for(int i=0;i<36;++i)
			    imp.get(a);
			 if(enter==false)
			 {
				 dspmatx.push_row(0);
				 rownum=dspmatx.getrnum()-1;
				 while(imp.get(a) && a!='\n')
					 dspmatx.push_to_row(rownum,a);
			 }
			 else if(enter==true)
			 {
				  while(imp.get(a) && a!='\n')
					  dspmatx.push_to_row(rownum,a);
				  rownum+=1;
			 }
		 }
		 else if(a==' ')
		 {
			 enter=true;
			 rownum=0;
			 for(int i=0;i<35;++i)
				 imp.get(a);
			 while(imp.get(a) && a!='\n')
			 {
				 identical.pushback(a);
			 }
		 }
	 }
	 dspmatx.writecolnum(0);                          //the "_column" member of "dspmatx" is zero unless a number is assigned
}
//Menu:[PFSP]-[_LogRLPmatx]
//Calculate log relative probability of MSA.(See Science 286, 295 (1999);for detail)
//Take 3 values."_seqsidmatx"aligned AAid seqs."_posi"position where mutation occur,-1 means no mutation."_AAid" AA id of the AA mutate to.
//Output 2 value."logrlpmatx" log relative p matrix."posiAAcnt" a 20*col darray holds number of AA at each position
void _LogRLPmatx(darray<int> & _seqsidmatx,darray<double> & logrlpmatx,darray<int> & posiAAcnt,const int & _posi=-1,const int & _AAid=-1)
{
	darray<int> tmpseqsidmatx;
	if(_posi<0)                                                                           //if no mutation occur,assign the original "_seqsidmatx"
		tmpseqsidmatx=_seqsidmatx;
	else if(_posi>=0)                                                                     //if mutate to "_AAid",locate the sequence and add to "tmpseqsidmatx"
	{
		for(int i=0;i<_seqsidmatx.getrnum();++i)
		{
			if(_seqsidmatx(i,_posi)==_AAid)
				tmpseqsidmatx.push_row(_seqsidmatx[i]);
		}
	}
	int rownum=tmpseqsidmatx.getrnum(),colnum=tmpseqsidmatx.getcnum(),positotalcnt;
	posiAAcnt.fast_resize(20+1,colnum+1,0);                                              //20 rows hold number count for 20 kinds of AA and the 21th row holds the sum value
	logrlpmatx.fast_resize(20,colnum);                                                   //log relative p matrix 
	for(int j=0;j<colnum;++j)                                                            //calculate AA count for each row and store in "posiAAcnt"
	{
		positotalcnt=0;
		for(int i=0;i<rownum;++i)
		{
			if(tmpseqsidmatx(i,j)>=0)                                                     
			{
				posiAAcnt(tmpseqsidmatx(i,j),j)+=1;
                positotalcnt+=1;
			}
		}
        posiAAcnt(20,j)=positotalcnt;
		for(int k=0;k<21;++k)
			posiAAcnt(k,colnum)+=posiAAcnt(k,j);
	}
	double p;
	bool found;
	for(int j=0;j<colnum;++j)                                                                            //calculate log relative p   
	{
		for(int i=0;i<20;++i)
		{
			found=false;
			for(int k=0;k<j;++k)
			{
				if(posiAAcnt(i,j)==posiAAcnt(i,k) && posiAAcnt(20,j)==posiAAcnt(20,k))
				{
					logrlpmatx(i,j)=logrlpmatx(i,k);
					found=true;
					break;
				}
			}
	        if(found==false)
			{
			   p=(double)posiAAcnt(i,colnum)/posiAAcnt(20,colnum);                                        //relative p
               if(p==0)
                  logrlpmatx(i,j)=0;                                                                                                   
			   else
			      logrlpmatx(i,j)=LogBinoP(rownum,posiAAcnt(i,j),p)-log(p);                               //equal to log(BinoP(...)/p)
			}
		}
	}
	return;
}
//Menu:[PFSP]-[_CalDetaG]
//Calculate coupling energy for each site
//Take 1 value."_logrlpmatx"log relative p matrix
//Output 1 value."DGarr"detaG (coupling energy) array
void _CalDetaG(darray<double> & _logrlpmatx,sarray<double> & DGarr)
{
	double temp;
	DGarr.resize(_logrlpmatx.getcnum());                               
	for(int j=0;j<_logrlpmatx.getcnum();++j)
	{
		temp=0;
		for(int i=0;i<20;++i)
			temp+=pow(_logrlpmatx(i,j),2.0);
		DGarr[j]=pow(temp,0.5);
	}
}
//Menu:[PFSP]
//Calculate DetaG and Deta_DetaG and output detail result
//Take 4 values."_seqsidmatx" aligned AAid seqs."_displaymatx" aligned AA seqs."_identical" sarray<char> store identical AA information."_seqnumthld" sequence number threshold
void PFSP(darray<int> & _seqsidmatx,darray<char> &_dspmatx,sarray<char> &_identical,const int & _seqnumthld)
{  
	clock_t start,end;
    start=clock();
    darray<int> posiAAcntmatx,tempmatx;                                                  //matrix holds AA count at each position
	darray<double> logrlpmatxtotal,logrlpmatxsub;                                        //log relative p matrix total/sub
	sarray<double> DGarr,tmpDGarr,D_DGtotal;                                             //DetaG array and Deta_DetaG array
	cout<<"Estimate conservative sites..."<<endl;
	_LogRLPmatx(_seqsidmatx,logrlpmatxtotal,posiAAcntmatx);                              //get log relative p matrix for the original matrix
	_CalDetaG(logrlpmatxtotal,DGarr);                                                    //compute detaG of each position and store in "DGarr"
	tmpDGarr=DGarr;                                                                      //assign value to temp var "tmpDGarr",which will be unrecoverably sorted
	QuickSort(tmpDGarr,0,tmpDGarr.size()-1);                                             //get detaG 95th/90th quartile
	double thldDG95=tmpDGarr[(int)((double)tmpDGarr.size()*95.0/100.0-1)];
	double thldDG90=tmpDGarr[(int)((double)tmpDGarr.size()*90.0/100.0-1)];
	funsitedata ProFunSite;                                                              //struct hold all result of coupling energy analysis
	int tmp=0;
    sarray<char> tmpAA;
	sarray<double> tmpD_DGarr;
	darray<double> tmpD_DGmatx;
	cout<<"Compute coupling energy...";
	for(int j=0;j<DGarr.size();++j)
	{                                                             
		PctMarker(j,DGarr.size()-1,2);
		for(int i=0;i<20;++i)
		{
           if(posiAAcntmatx(i,j)>=_seqnumthld && (posiAAcntmatx(20,j)-posiAAcntmatx(i,j)>=_seqnumthld))  //at list "_seqnumthld" seqs are required to make the statistics reasonable
		   {
			   tmpAA.pushback(idaa(i));                                                 //mutatable AA at each position
			   _LogRLPmatx(_seqsidmatx,logrlpmatxsub,tempmatx,j,i);
               _CalDetaG(logrlpmatxtotal-logrlpmatxsub,tmpD_DGarr);                     //coupling energy of other position when position "j" mutated to "i"
			   tmpD_DGmatx.push_row(tmpD_DGarr);
		   }
		}
		if(tmpAA.size()>0)                                                              //store all information in ProFunSite
		{
			ProFunSite._funsitearr.pushback(j);
			ProFunSite._AAvec.push_back(tmpAA);
			ProFunSite._D_DGvec.push_back(tmpD_DGmatx);
			tmp+=tmpAA.size();
			tmpAA.clear();
			tmpD_DGmatx.clear();
		}
	}
	cout<<endl<<"Compute coupling energy threshold...";

	int I=0,J=0;
	double tmpD_DG=0;
	sarray<double> tmparr1,tmparr2(tmp);                                   //"tmparr1" stores average D_DG for struct "ProFunSite"."tmparr2" also stores average D_DG,it is used for QuickSort
	D_DGtotal.resize(tmp*_seqsidmatx.getcnum());                           //pre compute the size of all D_DG data
	for(int i=0;i<ProFunSite._D_DGvec.size();++i)                                          
	{
		PctMarker(i,ProFunSite._D_DGvec.size()-1,2);
		tmparr1.resize(ProFunSite._D_DGvec[i].getrnum());                  //resize the size of "temparr1" to the number of mutatable AA at each position
		for(int j=0;j<ProFunSite._D_DGvec[i].getrnum();++j)
		{
			for(int k=0;k<ProFunSite._D_DGvec[i].getcnum();++k)
			{
				D_DGtotal[I]=ProFunSite._D_DGvec[i](j,k);                  //"D_DGtotal" stores all deta_deta G
				if(k!=ProFunSite._funsitearr[i])
				   tmpD_DG+=ProFunSite._D_DGvec[i](j,k);
				I+=1;
			}
            tmparr2[J]=tmparr1[j]=tmpD_DG/(DGarr.size()-1);                //store average D_DG into "tmparr1" and "temparr2"
			tmpD_DG=0;
			J+=1;
		}
		ProFunSite._avgD_DGvec.push_back(tmparr1);
	}
	cout<<endl<<"Output result..."<<endl;
	QuickSort(D_DGtotal,0,D_DGtotal.size()-1);
	double thldD_DG99=D_DGtotal[(int)((double)D_DGtotal.size()*99.0/100.0-1)];
	double thldD_DG95=D_DGtotal[(int)((double)D_DGtotal.size()*95.0/100.0-1)];
	double thldD_DG90=D_DGtotal[(int)((double)D_DGtotal.size()*90.0/100.0-1)];
    
	QuickSort(tmparr2,0,tmparr2.size()-1);
	double thldavgD_DG99=tmparr2[(int)((double)tmparr2.size()*99.0/100.0-1)];
	double thldavgD_DG95=tmparr2[(int)((double)tmparr2.size()*95.0/100.0-1)];
    double thldavgD_DG90=tmparr2[(int)((double)tmparr2.size()*90.0/100.0-1)];
	end=clock();
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////   write log
	char time[10],date[10];
	_strtime(time);_strdate(date);
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<" Process time "<<end-start<<" minisec"<<endl;
	outlog<<"Evaluate protein functional sites from "<<_seqsidmatx.getrnum()<<" sequences"<<endl;
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
	int width=(int)(log((double)_dspmatx.getrnum())/log(10.0))+2;                                                 //dynamic ID marker length.control output format
	sarray<char> BLANK(width,' '),conserve(DGarr.size());
	for(int i=0;i<DGarr.size();++i)
	{
		if(DGarr[i]<thldDG90)
			conserve[i]='-';
		else if(DGarr[i]>=thldDG95)
			conserve[i]='h';
		else
			conserve[i]='m';
	}
	MKDIR("C:\\PFSP\\FUNCSITE");
	string path="C:\\PFSP\\FUNCSITE\\PFSP_"+inttostr(_dspmatx.getrnum())+"_"+inttostr(_seqnumthld)+".txt";
	ofstream outPFS(path.c_str());
	outPFS<<date<<" "<<time<<endl;
	outPFS<<"Protein Functional Sites Prediction--PFSP:(C++) Copyright 2009,Ruan Xiaoyang"<<endl;
	outPFS<<"Information drawn from "<<_seqsidmatx.getrnum()<<" sequences"<<endl<<endl;

    outPFS<<"Conservative site estimation (DetaG)"<<endl;
	outPFS<<"Conservation  '*':identical  '-':low(<90% percentile)  'm':moderate(90% percentile DetaG:"<<thldDG90<<"KT)  'h':high(95% percentile DetaG:"<<thldDG95<<"KT)  DetaG(min):"<<tmpDGarr[0]<<"KT)  DetaG(max):"<<tmpDGarr[tmpDGarr.size()-1]<<"KT)"<<endl;
	outPFS<<BLANK<<" ";LandMarker(DGarr.size(),outPFS);outPFS<<endl;
    outPFS<<BLANK<<" ";MileStone(DGarr.size(),outPFS);outPFS<<endl;
    outPFS<<BLANK<<" "<<conserve<<endl;
	outPFS<<BLANK<<" "<<_identical<<endl;
	for(int i=0;i<_dspmatx.getrnum();++i)
	{
		if(i!=0 && i!=(_dspmatx.getrnum()-1) && i%50==0)
		{
			outPFS<<BLANK<<" ";LandMarker(DGarr.size(),outPFS);outPFS<<endl;
            outPFS<<BLANK<<" ";MileStone(DGarr.size(),outPFS);outPFS<<endl;
            outPFS<<BLANK<<" "<<conserve<<endl;
		}
		outPFS.setf(ios::right);outPFS.width(width);
		outPFS<<i<<" "<<_dspmatx[i]<<" "<<i<<endl;
	}
    outPFS<<BLANK<<" "<<_identical<<endl;
	outPFS<<BLANK<<" "<<conserve<<endl;
	outPFS<<BLANK<<" ";MileStone(DGarr.size(),outPFS);outPFS<<endl;
	outPFS<<BLANK<<" ";LandMarker(DGarr.size(),outPFS);outPFS<<endl<<endl;
    
	double avgD_DGtmp;
    outPFS<<"Average coupling energy(average Deta_DetaG)"<<endl
		<<"' 'data not available  '-':LOW(<90% percentile)  'M':moderate(90% percentile AvgDeta_DetaG:"<<thldavgD_DG90<<"KT)  'H':high(95% percentile AvgDeta_DetaG:"<<thldavgD_DG95<<"KT)  'E':extreme(99% percentile AvgDeta_DetaG:"<<thldavgD_DG99<<"KT)"<<endl;
	outPFS<<BLANK<<" ";LandMarker(DGarr.size(),outPFS);outPFS<<endl;
	outPFS<<BLANK<<" ";MileStone(DGarr.size(),outPFS);outPFS<<endl;
    outPFS<<BLANK<<" ";
	for(int i=0;i<ProFunSite._avgD_DGvec.size();++i)
	{
		for(int j=0;j<(ProFunSite._funsitearr[i]-(i==0?-1:ProFunSite._funsitearr[i-1])-1);++j)
			outPFS<<" ";
		avgD_DGtmp=ProFunSite._avgD_DGvec[i].smax();
        if(avgD_DGtmp<thldavgD_DG90)
			outPFS<<"-";
		else if(avgD_DGtmp<thldavgD_DG95)
			outPFS<<"M";
		else if(avgD_DGtmp<thldavgD_DG99)
			outPFS<<"H";
		else
			outPFS<<"E";
	}
	outPFS<<endl<<endl;
	
	outPFS<<"Statistical coupling energy estimation (Deta_DetaG)"<<endl;
	outPFS<<"Coupling energy  '-':low(<90% percentile)  'm':moderate(90% percentile Deta_DetaG:"<<thldD_DG90<<"KT)  'h':high(95% percentile Deta_DetaG:"<<thldD_DG95<<"KT)  'e':extreme(99% percentile Deta_DetaG:"<<thldD_DG99<<"KT)"<<endl;
	outPFS<<"Deta_DetaG(min):"<<D_DGtotal[0]<<"KT)  Deta_DetaG(max):"<<D_DGtotal[D_DGtotal.size()-1]<<"KT)"<<endl;
	double D_DGtmp;
	for(int n=0;n<ProFunSite._D_DGvec.size();++n)
	{
		outPFS<<"Site "<<ProFunSite._funsitearr[n]<<"  ";
		for(int m=0;m<ProFunSite._AAvec[n].size();++m)
		{
			outPFS<<"->"<<ProFunSite._AAvec[n][m]<<"(Avg "<<ProFunSite._avgD_DGvec[n][m]<<"KT ";
			 if(ProFunSite._avgD_DGvec[n][m]<thldavgD_DG90)
			outPFS<<"'LOW')  ";
		     else if(ProFunSite._avgD_DGvec[n][m]<thldavgD_DG95)
			outPFS<<"'M')  ";
		     else if(ProFunSite._avgD_DGvec[n][m]<thldavgD_DG99)
			outPFS<<"'H')  ";
		     else
			outPFS<<"'E')  ";
		}
		outPFS<<endl;
	    outPFS<<BLANK<<" ";LandMarker(DGarr.size(),outPFS);outPFS<<endl;
	    outPFS<<BLANK<<" ";MileStone(DGarr.size(),outPFS);outPFS<<endl;
		for(int i=0;i<ProFunSite._D_DGvec[n].getrnum();++i)
		{
		    outPFS<<BLANK<<" ";
			for(int j=0;j<ProFunSite._D_DGvec[n].getcnum();++j)
			{
                D_DGtmp=ProFunSite._D_DGvec[n](i,j);
				if(j==ProFunSite._funsitearr[n])
				{
					outPFS<<ProFunSite._AAvec[n][i];
					continue;
				}
				if(D_DGtmp<thldD_DG90)
					outPFS<<"-";
				else if(D_DGtmp>=thldD_DG99)
					outPFS<<"e";
				else if(D_DGtmp>=thldD_DG95)
					outPFS<<"h";
				else
					outPFS<<"m";
			}
			outPFS<<endl;
		}
	}
	outPFS.close();
	cout<<"PFSP Result has been saved to "<<path<<endl;
	TMPPATH=path;
	DSPFILE(TMPPATH);
}
//Compute funtional site from MSA result
//Take 2 value."_path"the path where MSA result was saved (Note:the 1st row of MSA store row and column number,the 2nd row store identical information '*'identical ' 'non-identical)."_seqnumthld" sequence number threshold
//Return 0 if MSA has less than 20 seqs.Return 1 if successful.
int PFSP(char * _path,const int & _seqnumthld)
{
	darray<char> dspmatx;
	sarray<char> identical;
	cout<<"Specify [MSA result format]\n"
		<<"ep:P represent [PFSP MSA result]\n"
		<<"   C represent [ClustalW MSA result]\n"
		<<"Input:";
	char MSAformat;
    cin>>MSAformat;
	cout<<"Import file..."<<endl;
	if(MSAformat=='p'|| MSAformat=='P')
	    _ImpMSA(_path,dspmatx,identical);                                     //import file to identical array and display matrix
	else if(MSAformat=='c'|| MSAformat=='C')
		_ImpMSAClustalW(_path,dspmatx,identical);
	if(dspmatx.getrnum()<20)                                              //compute coupling energy for MSA with less than 20 seqs is somewhat useless.the core idea is theoretically every AA has a chance to appear in 20 random selections
	   return 0;
	darray<int> seqsidmatx(dspmatx.getrnum(),dspmatx.getcnum());          //transform to AAid seqs
	for(int i=0;i<dspmatx.getrnum();++i)
	{
		for(int j=0;j<dspmatx.getcnum();++j)
			seqsidmatx(i,j)=aaid(dspmatx(i,j));
	}
	PFSP(seqsidmatx,dspmatx,identical,_seqnumthld);
    return 1;
}



