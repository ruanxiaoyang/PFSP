//////////////////////////////////////////////////////////////
//                        RBLAST                            //
//           Copyright Belongs to Ruan Xiaoyang             //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
//This module performs glocal(global-local) alignment between two amino acid sequences (not for nucleotide sequence)
//Word matching rule was first applied to find local alginment,then Needleman global alignment rule was applied to align the remaining part
//For each run,a log file was saved at c:\ProFunSite.log indicating the alignment information


struct wordlib                                                     //structure hold sequence word library
{
	darray<int> _word;                                            
	darray<int> _location;
	int _libsize;
};

//inline void T_Sthld();
//inline void _CRlib();
//inline void _Extend();
//inline void _Converge();
//inline void _Terminal();
//inline darray<int> _Blankmark();
//inline void _Fillblank();
//inline void _Maxscwp();
//template<class type>
//inline type _Getscore();
//menu;[RBLAST]-[BLAST]-[T_Sthld]
//                     -[_CRlib]
//                     -[_Extend]
//                     -[_Converge]-[_Terminal]
//             -[_Blankmark]
//             -[_Fillblank]-[_Maxscwp]

//"_wlen" word length. _wlen=n allows a maximum of n-1 consecutive negative score mismatch.larger genetic distance should use higher _wlen to allow for longer mismatch
//"_scorematx" score matrix
//"_qery" quiry sequence
//"_bseq" base sequence
//"_scmatxid" score matrix id
//"_gappnt" gap penalty
//"_cutoffper" cutoff percentile.it defined the t and s score threshold above which a hit and a maximum alignment will be admitted 
//"_DorB" determined the AA frequency database (Dayhoff or Blosum) used to generate random AA sequence
//compute "score"
//compute standard 0 or 1 alignment matrix "_matrix"
template<class type>                                                 
void RBLAST(const int & _wlen,darray<type> & _scorematx,sarray<int> & _qeryid,sarray<int> & _bseqid,const int & _scmatxid,darray<int> & _matrix,type & score,const int &_gappnt=-3,const int & _cutoffper=95,const char &_DorB='D')
{
	clock_t start,end;
	start=clock();
	darray<int> blank;
	BLAST(_wlen,_scorematx,_qeryid,_bseqid,_scmatxid,_matrix,_cutoffper,_DorB);
	_Blankmark(_matrix,_qeryid.size(),_bseqid.size(),blank);                                      //find blank area not covered by word matching rule
	_Fillblank(blank,_qeryid,_bseqid,_matrix);                                                    //fill in blank area by Needleman global alignment rule
    score=_Getscore(_scorematx,_matrix,_qeryid.size(),_bseqid.size(),_qeryid,_bseqid,_gappnt);
	end=clock();
	////////////////////////////////////////////////////////////////////////////////////////////write log
	char time[10],date[10];
	_strtime(time);_strdate(date);
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<" Process time "<<end-start<<" minisec"<<endl;
	outlog<<"Glocal alignment of "<<_qeryid.size()<<" AA seq with "<<_bseqid.size()<<" AA seq."<<endl;
	outlog<<"Word length:"<<_wlen<<" Score matrix:"<<(_DorB=='D' ? "PAM": "Blosum")<<_scmatxid<<" Cutoff percentile:"<<_cutoffper<<endl<<endl; 

}
template <class type>
void RBLAST(proformat & _seqsandinf,vector<darray<double>> &_scoredb,darray<type> & _scorematx,const int & _scmatxid,const char & _fixscmatx,const char & _DorB,const char & _EXorAP,const char & _subtype,const int & _gappnt=-3,const int & _gma=1,const int & _wlen=3,const int & _step=10,const int & _cutoffper=95,const int & _acculvl=4)
{
	MKDIR("C:\\PFSP\\RBLAST");
	string cmnpath="c:\\PFSP\\RBLAST\\RBLAST_"+inttostr(_seqsandinf._seqnum)+"_"+string(&_fixscmatx,sizeof(char))+"_"+string(&_DorB,sizeof(char));
	string path;
	if(_fixscmatx=='Y'||_fixscmatx=='y')
	   path=cmnpath+"_"+inttostr(_scmatxid)+"_"+inttostr(_gappnt)+".txt";
	else
	   path=cmnpath+"_"+string(&_EXorAP,sizeof(char))+"_"+string(&_subtype,sizeof(char))+"_"+inttostr(_gappnt)+".txt";
	ofstream outRBL(path.c_str());
	char time[10],date[10];
	_strtime(time);_strdate(date);
    outRBL<<date<<" "<<time<<endl;
	outRBL<<"Protein Functional Sites Prediction--RBLAST:(C++) Copyright 2009,Ruan Xiaoyang"<<endl;
	outRBL<<"Pairwise align "<<_seqsandinf._seqnum<<" sequences   Word length:"<<_wlen<<"  Gap penalty:"<<_gappnt<<"  ";
	if(_fixscmatx=='Y' || _fixscmatx=='y')
		outRBL<<"Scoring matrix:"<<((_DorB=='D'|| _DorB=='d')?"PAM":"BLOSUM")<<_scmatxid<<endl<<endl;
	else
	{
		if(_EXorAP=='A'||_EXorAP=='a')
		{
			outRBL<<"Score database:DayhoffPAM1-400  Approximate distance:"<<_subtype;
			if(_subtype=='G' || _subtype=='g')
				outRBL<<"  Gamma:"<<_gma;
			outRBL<<endl<<endl;
		}
		else if(_EXorAP=='E'||_EXorAP=='e')
			outRBL<<"Score database:DayhoffPAM1-400  Bayesian exact distance:"<<_subtype<<"  Bayesian step:"<<_step<<endl<<endl;
	}
	type score;double dist,var,smlt;
	darray<int> alnmatx;
	darray<type> scorematx;
	int width=(int)(log((double)_seqsandinf._seqnum)/log(10.0))+2,scmatxid;
	cout<<"RBLAST pairwise alignment in progress...";
	for(int i=0;i<_seqsandinf._seqnum;++i)
	{
		PctMarker(i,_seqsandinf._seqnum-1,2);
        for(int j=i+1;j<_seqsandinf._seqnum;++j)
		{
			if(_seqsandinf._format=='C')                               //when using custom format,only output _gi(stored all information inputed by user
			{
				    outRBL<<"ID";
					outRBL.width(width);
					outRBL<<j<<" "<<_seqsandinf._infvec[j]._gi<<endl;
					outRBL<<"ID";
					outRBL.width(width);
					outRBL<<i<<" "<<_seqsandinf._infvec[i]._gi<<endl;
			}
			else if(_seqsandinf._format=='F')
			{
					outRBL<<"ID";
					outRBL.width(width);
					outRBL<<j<<" Gi_number:"<<_seqsandinf._infvec[j]._ginum<<" Accession:";
					outRBL<<_seqsandinf._infvec[j]._accession<<" Locus:"<<_seqsandinf._infvec[j]._locus<<endl;
					outRBL<<"ID";
					outRBL.width(width);
					outRBL<<i<<" Gi_number:"<<_seqsandinf._infvec[i]._ginum<<" Accession:";
					outRBL<<_seqsandinf._infvec[i]._accession<<" Locus:"<<_seqsandinf._infvec[i]._locus<<endl;
			}
			if(_fixscmatx=='Y' || _fixscmatx=='y')
			{
				scorematx=_scorematx;
                scmatxid=_scmatxid;
			}
			else
			{
				if(_EXorAP=='A'||_EXorAP=='a')
				{
				   darray<int> BKMATX;
				   AprxDist(_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],dist,var,smlt,_subtype,BKMATX,_gma);
				}
				else if(_EXorAP=='E'||_EXorAP=='e')
				{
					BayesianDist(_subtype,_scoredb,_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],dist,var,_wlen,_gappnt,_step,_acculvl);
					if(_subtype=='B' || _subtype=='b')
						outRBL<<"(Local)";
					else
						outRBL<<"(Global)";
				}
			    scorematx=_scoredb[(int)((dist*100)-1)];
				scmatxid=(int)(dist*100);
				outRBL<<"Estimated distance:PAM"<<(int)(dist*100)<<"  Variance:"<<var*100<<"  ";
				if(_EXorAP=='A'||_EXorAP=='a')
					outRBL<<"Similarity:"<<smlt*100<<"%   ";

			}
			RBLAST(_wlen,scorematx,_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],scmatxid,alnmatx,score,_gappnt,_cutoffper,_DorB);
			outRBL<<"Score(10*log[base10]):"<<score<<endl;
			Dsptwo(_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],scorematx,alnmatx,outRBL);
			outRBL<<endl;
		}
	}
	cout<<endl;
	outRBL.close();
	cout<<"Result has been saved to "<<path<<endl;
	writeTStofile();
	TMPPATH=path;
	DSPFILE(TMPPATH);
}

struct CRlibstr
{
	int _wlen;
	sarray<int> * _seqid;
	wordlib * _seqlib;
};

CRITICAL_SECTION WRITETS;       //This var control access to global var "TSSTRVEC",used in BLAST
HANDLE LIBSTRREADY;             //This var is used by CreateEvent function in BLAST to control access to CRlibstr
template<class type>                                                 
void BLAST(const int & _wlen,darray<type> & _scorematx,sarray<int> & _qeryid,sarray<int> & _bseqid,const int & _scmatxid,darray<int> & _matrix,const int & _cutoffper=95,const char &_DorB='D')
{
	bool recfound=false;
	type tscore,tthld,sthld;
	int rownum=_qeryid.size(),colnum=_bseqid.size(),seqlen=(rownum>colnum ? colnum:rownum),err;
	err=((int)((double)seqlen/200.0)+1)*5;               //when sequence length is large,there is generally no need to be very accurate
	for(int i=0;i<TSSTRVEC.size();++i)                   //look into existing records for match to save time
	{
		if(_DorB==TSSTRVEC[i]._DorB && _scmatxid==TSSTRVEC[i]._scmatxid && _wlen==TSSTRVEC[i]._wlen && TSSTRVEC[i]._seqlen>=seqlen-err && TSSTRVEC[i]._seqlen<=seqlen+err && TSSTRVEC[i]._cutoffper==_cutoffper) //sequence length is acceptable within a short range
		{
			tthld=TSSTRVEC[i]._tthld;
			sthld=TSSTRVEC[i]._sthld;
			recfound=true;
			break;
		}
	}
	if(recfound==false)                                   //if no record is found,we have to engage T_Sthld function
	{
		T_Sthld(seqlen,_wlen,_cutoffper,_scorematx,tthld,sthld,_DorB);                              //if no record found,calculate t and s score threshold and store results
	    _TSSTR temp;
		temp._DorB=_DorB;
		temp._scmatxid=_scmatxid;
		temp._wlen=_wlen;
		temp._seqlen=seqlen;
		temp._cutoffper=_cutoffper;
		temp._tthld=tthld;
		temp._sthld=sthld;
		InitializeCriticalSection(&WRITETS);               //this critical section controls simultanous access to global var "TSSTRVEC"
		EnterCriticalSection(&WRITETS);
		TSSTRVEC.push_back(temp);
		LeaveCriticalSection(&WRITETS);
	}
	darray<type> tempmatx(rownum,colnum,(type)0);                                                      //its necessary to initialize to 0 (for _Converge function)
	wordlib qerylib,bseqlib;
	if(BYSINPROGRESS==false)               //this part is used when [BysDistMatx],which will occupy all CPU resource,is not in progress 
	{
		LIBSTRREADY=CreateEvent(NULL,FALSE,FALSE,NULL);
		CRlibstr CRlibdata;
		CRlibdata._wlen=_wlen;
		CRlibdata._seqid=&_qeryid;
		CRlibdata._seqlib=&qerylib;
		HANDLE handle[2];
		handle[0]=(HANDLE)_beginthreadex(NULL,0,&_CRlib,&CRlibdata,NULL,NULL);                         //create word library for quiry and base sequences
		WaitForSingleObject(LIBSTRREADY,INFINITE);
		CRlibdata._seqid=&_bseqid;
		CRlibdata._seqlib=&bseqlib;
		handle[1]=(HANDLE)_beginthreadex(NULL,0,&_CRlib,&CRlibdata,NULL,NULL);
		WaitForMultipleObjects(2,handle,TRUE,INFINITE);
	}
	else if(BYSINPROGRESS==true)          //this part is used to avoid creating additional thread when [BysDistMatx] is triggered
	{
		_CRlib(_wlen,_qeryid,qerylib);
		_CRlib(_wlen,_bseqid,bseqlib);
	}
	for(int i=0;i<qerylib._libsize;++i)
	{
		for(int j=0;j<bseqlib._libsize;++j)
		{
			tscore=(type)0;
			for(int w=0;w<_wlen;++w)
				tscore+=_scorematx(qerylib._word(i,w),bseqlib._word(j,w));                             //calculate score of hit
			if(tscore>=tthld)                                                                          //if "tscore" reached threshold then extend hit
			{ 
				for(int k=0;k<qerylib._location.getcnum(i);++k)
				{
					for(int m=0;m<bseqlib._location.getcnum(j);++m)
						_Extend(_wlen,sthld,qerylib._location(i,k),bseqlib._location(j,m),_scorematx,_qeryid,_bseqid,tempmatx);
				}
			}
		}
	}
    _Converge(tempmatx,rownum,colnum);                                               //converge way point to get unique alignment pathway
	_matrix.fast_resize(rownum,colnum);
	for(int i=0;i<rownum;++i)                                                      //retore alignment matrix to standard 0 or 1 form
	{
		for(int j=0;j<colnum;++j)
		{
			if(tempmatx(i,j)==0 || tempmatx(i,j)==INT_MIN)
				_matrix(i,j)=0;
			else
				_matrix(i,j)=1;
		}
	}
}

//Menu:[RBLAST]-[BLAST]-[T_Sthld]
//Calculate t and s score distribution for two random sequences.
//Take 5 values."_seqlen" sequence length of the shorter sequence."_wlen" word length."_cutoffper" cutoff percentile."_scmatx" score matrix used."_DorB" Dayhoff or Blosum AA frequency database 
//Compute 2 values."tthld" t score threshold."sthld" s score threshold.
template <class type>  
void T_Sthld(const int & _seqlen,const int & _wlen,const int & _cutoffper,darray<type> & _scmatx,type & tthld,type & sthld,char _DorB='D')
{
	int num=5000,mark=0;                                              //"num" number of replication
	type tscore,sscore,tempscore;
	sarray<int> seqa,seqb;
	sarray<type> scorevec(num,(type)0),temparr;
	sarray<int> accumdb;
	if(_DorB=='D')
	   accumdb=DayAccAAFreq;                                          //accumulate AA frequency database.used for generate random seq.
	else
       accumdb=BloAccAAFreq;
    for(int n=0;n<num;++n)
	{
		tscore=(type)0;
		randAAidseq(_wlen,accumdb,seqa);                              //generate random AA seq with "_wlen" length
		randAAidseq(_wlen,accumdb,seqb);
		for(int w=0;w<_wlen;++w)
			tscore+=_scmatx(seqa[w],seqb[w]);
		scorevec[n]=tscore;
	}
	QuickSort(scorevec,0,num-1);                                      //sort "scorevec" from min to max
    tthld=scorevec[(int)((double)num*(double)_cutoffper/100.0)-1];    
	for(int n=0;n<num;++n)
	{
		randAAidseq(_seqlen,accumdb,seqa);
		randAAidseq(_seqlen,accumdb,seqb);
		sscore=(type)INT_MIN;                                         //initiate "sscore" to a very small value
		for(int i=0;i<_seqlen;++i)
		{
			if(_scmatx(seqa[i],seqb[i])>=0)                                             //only start pushback when meet a positive score
			{
				mark=0;
				temparr.pushback(_scmatx(seqa[i],seqb[i]));
			}
            else if(_scmatx(seqa[i],seqb[i])<0 && temparr.size()>0)                    
			{
				mark+=1;
			    temparr.pushback(_scmatx(seqa[i],seqb[i]));
			}
			if(temparr.size()>0 && (mark==_wlen || i==_seqlen-1))
			{
				tempscore=(type)0;
				temparr.resize(temparr.size()-mark);
				for(int j=0;j<temparr.size();++j)
					tempscore+=temparr[j];
		        if(tempscore>sscore)                                  //find the highest score in each pair of random alignment and store in "sscore"
					sscore=tempscore;
				temparr.clear();
			}
		}
		scorevec[n]=sscore;                                           
	}
	QuickSort(scorevec,0,num-1);
	sthld=scorevec[(int)((double)num*(double)_cutoffper/100.0)-1];
}

//Menu:[RBLAST]-[BLAST]-[_CRlib]
//Child thread creates word library for AA sequence (also capable with nucleotide)
//Take pointer to CRlibstr
//Fill in 1 "_seqlib" member in CRlibstr."seqlib" seuquence library structure, contain 3 members "seqlib._word" hold word information,"seqlib._location" hold location of word,"seqlib._libsize" hold size of "seqlib" 
unsigned __stdcall _CRlib(void * param)
{
	CRlibstr *CRlibdata=static_cast<CRlibstr*>(param);
	int k,wlen=CRlibdata->_wlen;
	sarray<int> * seqid=CRlibdata->_seqid;
	wordlib * seqlib=CRlibdata->_seqlib;
	SetEvent(LIBSTRREADY);
	sarray<int> temp(wlen,0);
	seqlib->_libsize=0;
	for(int i=0;i<seqid->size()-wlen+1;++i)
	{
		for(int w=0;w<wlen;++w)
			temp[w]=seqid->operator [](i+w);
		for(k=0;k<seqlib->_libsize;++k)
		{
			if(seqlib->_word[k]==temp)
				break;
		}
		if(k==seqlib->_libsize)                                      //if no previous record found,add a new line to _word,add a new location to _location
		{
			seqlib->_word.push_row(temp);
			seqlib->_location.push_row(1);
			seqlib->_location(seqlib->_libsize,0)=i;
		    seqlib->_libsize+=1;
		}
		else
			seqlib->_location.push_to_row(k,i);
	}
	return 0;
}
//Menu:[RBLAST]-[BLAST]-[_CRlib]
//Create word library for AA sequence (also capable with nucleotide)
//Take 2 values."_wlen" word length."_seq" AA sequence.
//Output 1 value."seqlib" seuquence library structure, contain 3 members "seqlib._word" hold word information,"seqlib._location" hold location of word,"seqlib._libsize" hold size of "seqlib" 
void _CRlib(const int & _wlen,sarray<int> & _seqid,wordlib & seqlib)
{
	int k;
	sarray<int> temp(_wlen,0);
	seqlib._libsize=0;
	for(int i=0;i<_seqid.size()-_wlen+1;++i)
	{
		for(int w=0;w<_wlen;++w)
			temp[w]=_seqid[i+w];
		for(k=0;k<seqlib._libsize;++k)
		{
			if(seqlib._word[k]==temp)
				break;
		}
		if(k==seqlib._libsize)                                      //if no previous record found,add a new line to _word,add a new location to _location
		{
			seqlib._word.push_row(temp);
			seqlib._location.push_row(1);
			seqlib._location(seqlib._libsize,0)=i;
		    seqlib._libsize+=1;
		}
		else
			seqlib._location.push_to_row(k,i);
	}
}

//Menu:[RBLAST]-[BLAST]-[_Extend]
//Extend hit to see if max score exceed threshold.
//Take 7 values."_wlen" word length."_sthld" s threshold."_qst"/"_bst" start position of quiry/base sequence."_scmatx" score matrix."_qery"/"_bseq" quiry/base sequence.
//Update "_matrix",alignment matrix
template<class type>                                                  
void _Extend(const int & _wlen,const double & _sthld,const int & _qst,const int & _bst,darray<type> & _scmatx,sarray<int> & _qeryid,sarray<int> & _bseqid,darray<type> & _matrix)
{
	type score;
	int checker=0,qery,bseq,tempqst=_qst,tempbst=_bst,qst,bst,right=0,left=0;          //"tempqst"/"tempbst" record the start position."left"/"right"record the left/right terminal where "_wlen" consecutive negative scores accur
	while(tempqst<_qeryid.size()-1 && tempbst<_bseqid.size()-1 && _qeryid[tempqst]!=_bseqid[tempbst])   //if score at start position <0, then move to right until score>=0
	{
		tempqst+=1;
		tempbst+=1;
	}
	if(_matrix(tempqst,tempbst)!=0)                                       //check existing record to see if this position has already been extended.if yes,then stop
		return;
	sarray<type> scorevec;                                            
	qst=tempqst;bst=tempbst;    
	while(qst>=0 && bst>=0)                                              //go backward until "_wlen" consecutive negative scores accur or reach the left terminal
	{
		qery=_qeryid[qst];
		bseq=_bseqid[bst];
		if(_scmatx(qery,bseq)>=0)
			checker=0;
		else
			checker+=1;
		left-=1;
		if(checker<=_wlen)
			scorevec.pushback(_scmatx(qery,bseq));
	    if(checker==_wlen || qst==0 || bst==0)
		{
			left+=checker;
			scorevec.resize(scorevec.size()-checker);                      //truncate negative score at the left end 
			break;
		}
		qst-=1;bst-=1;
	}
	left+=1;                                                       //left add 1 because the search started from tempqst and tempbst
	qst=tempqst+1;bst=tempbst+1;checker=0;                         //restore iterator to 1 pass start position
	while(qst<_qeryid.size() && bst<_bseqid.size())                        //go forward until "_wlen" consecutive negative scores accur or reach the right terminal
	{
		qery=_qeryid[qst];
		bseq=_bseqid[bst];
		if(_scmatx(qery,bseq)>=0)
			checker=0;
		else
			checker+=1;
	    right+=1;
		if(checker<=_wlen)
		    scorevec.pushback(_scmatx(qery,bseq));
		if(checker==_wlen || qst==_qeryid.size()-1 || bst==_bseqid.size()-1)
		{
			right-=checker;
			scorevec.resize(scorevec.size()-checker);                      //truncate negative score at the right end
			break;
		}
		qst+=1;bst+=1;
	}
	score=(type)0;
	for(int i=0;i<scorevec.size();++i)
        score+=scorevec[i];
	if(score>=_sthld)                                                //if alignment score is larger than threshold then update match matrix information
	{
        for(int i=left;i<=right;++i)
			_matrix(tempqst+i,tempbst+i)=score;
	}
}
//Menu:[RBLAST]-[BLAST]-[_Converge]
//Merge deviated/weak score alignment to better score alignment to obtain unique waypoint
//Take 3 values."_matrix" alignment matrix."_rownum"/"_colnum" row/colnum number of alignment matrix
//Update "_matrix"
template <class type>
void _Converge(darray<type> & _matrix,const int & _rownum,const int & _colnum)
{
	type comp1,comp2;
	int topI,topJ,botI,botJ,topi,topj,boti,botj;
	for(int i=0;i<_rownum;++i)                           
	{
		for(int j=_colnum-1;j>=0;--j)                                //iterate from right to left                              
		{
			if(_matrix(i,j)==0)                                      //if start point "_matrix(i,j)" equal to 0,move to next column
				continue;
			else if(_matrix(i,j)==INT_MIN)                           //if start point is row stop marker INT_MIN,move to next row
				break;
			else
			{
				comp1=_matrix(i,j);
				_Terminal(_matrix,i,j,topi,topj,boti,botj);
	            for(int ii=i;ii<_rownum;++ii)
	            {
		            for(int jj=j-1;jj>=0;--jj)
		            {
						if((comp2=_matrix(ii,jj))==INT_MIN)          //if row stop marker appear,check next row.speed up converge process   
							break;
			            if(comp2!=0)                                 //erase the part of weaker alignment that located in "dead zone" of stronger alignment
			            {
				            _Terminal(_matrix,ii,jj,topI,topJ,botI,botJ);
				            if(comp2>comp1 || (comp2==comp1 && abs(ii-jj)<=abs(i-j)))   //if the newly found segment is stronger than the target
				            {
				                for(int m=0;m<=boti-topi;++m)
				                {
					                 if((topi+m)<=botI && (topj+m)>=topJ)
						                 _matrix(topi+m,topj+m)=0;
				                }
				            }
				            if(comp2<comp1 ||  (comp2==comp1 && abs(ii-jj)>abs(i-j)))
			            	{
			               	   for(int m=0;m<=botI-topI;++m)
				               {
					            	if((topI+m)>=topi && (topJ+m)<=botj)
						        	    _matrix(topI+m,topJ+m)=0;
				               }
				             }
				             if(_matrix(i,j)==0)                    //only when start point has been erased can move to next column
	                			  goto nextj;
						}
					}
				}
				for(int ii=i;ii<_rownum;++ii)                       //if start point was preserved (!=0) after convergence,mark lower left area of start point as useless to gain speed boost
	            {
		            for(int jj=j-1;jj>=0;--jj)
					{
						if(_matrix(ii,jj)==INT_MIN)
							break;
						else
							_matrix(ii,jj)=INT_MIN;
					}
				}
			}
nextj:;
		}
	}
	                                                                //coordinate marked with INT_MIN will be restored to 0 in [RBLAST]
}
//Menu:[RBLAST]-[BLAST]-[_Converge]-[_Terminal]
//Calculate diagonal terminal coordinate.A terminal is reached when "_matrix(_i,_j)" evaluate to 0
//Take 3 values."_matrix" alignment matrix of quiry and base sequence."_i"/"_j" row/column id of start point.
//Compute 4 values."topi"/"topj" top terminal row/column id."boti"/"botj" bottom terminal row/column id. 
template<class type>
void _Terminal(darray<type> & _matrix,const int & _i,const int & _j,int & topi,int & topj,int & boti,int & botj)
{
	int i=_i-1,j=_j-1;
	while(i>=0 && j>=0 && _matrix(i,j)>0)                                                //go backward to find top terminal
	{
			i-=1;j-=1;
	}
	topi=i+1;topj=j+1;
	i=_i+1;j=_j+1;
	while(i<_matrix.getrnum() && j<_matrix.getcnum() && _matrix(i,j)>0)                //go forward to find bottom terminal
	{
			i+=1;j+=1;
	}
    boti=i-1,botj=j-1;
}
//Menu:[RBLAST]-[_Blankmark]
//Find out blank areas not covered by word extending strategy
//Take 3 values."alnmatx" alignment matrix (standard 0 or 1 form)."_rownum"/"_colnum" row/column number of "alnmatx"
//Teturn 4-column darry<int> containing coordinates marking the blank area.column 0,1,2,3 are top i,top j,bottom i,bottom j respectively
void _Blankmark(darray<int> & alnmatx,const int & _rownum,const int & _colnum,darray<int> & blank)
{
	bool open=false;
	int i=0,j=0,ii,jj,lr;                                //"lr" last row id
    while(i<_rownum && j<_colnum)
	{
	    for(jj=j;jj<_colnum;++jj)                        //scan the row for 1
		{
			if(alnmatx(i,jj)!=0)                         //if record found,record the column id and break
			{
				j=jj;
				break;
			}
		}
		if(jj==_colnum)                                  //only start scanning the column when no record found in the row
		{
			for(ii=i;ii<_rownum;++ii)
			{
				if(alnmatx(ii,j)!=0)                     //if record found,record the row id and break
				{
					i=ii;
					break;
				}
			}
		}
		if(jj==_colnum && ii==_rownum && open==false)    //when no record found in row and column and blank has not been openned yet,open a blank
		{
			open=true;
			blank.push_row(4);
			lr=blank.getrnum()-1;
			blank(lr,0)=i,blank(lr,1)=j;
		}
		else if((jj!=_colnum || ii!=_rownum) && open==true) //if a blank is open,and a record was found,close the blank and store the end coordinate
		{
				blank(lr,2)=i-1;
				blank(lr,3)=j-1;
				open=false;
		}
		i+=1;j+=1;
	}
	if(open==true)                                          //if the lower right terminal of "alnmatx" is 0,filled in the coordinate 
	{
		blank(lr,2)=_rownum-1;blank(lr,3)=_colnum-1;
	}
}

//Menu:[RBLAST]-[_Fillblank]-[_Maxscwp]
//Find score waypoint for blank area after [_Converge]
//Take 2 values."_BKMATX" Needleman global alignment matrix."wp" waypoint darray
//Update 1 value."wp" waypoint darray:column 0 store row id,column 1 store column id.
void _Maxscwp(darray<int> & _BKMATX,darray<int> & wp)
{
	int i=0,j=0,r,s,d,start,row=_BKMATX.getrnum()-1,col=_BKMATX.getcnum()-1;
	wp.fast_resize(0,2);
	bool ropen=false,dopen=false;
    while(i<row && j<col)
	{
		start=_BKMATX(i,j);                                                //scoer of the current cell
		r=_BKMATX(i,j+1);                                                  //score of the right cell
		d=_BKMATX(i+1,j);                                                  //score of the down cell
		s=_BKMATX(i+1,j+1);                                                //scoer of the diagonal cell
		     //when no gap open occur                   //when right gap open                    //when down gap open     
        if(((r!=start && d!=start) || (s>=r && s>=d))||(ropen==true && (r!=start || s>=r))||(dopen==true && (d!=start || s>=d)))
		{  
			ropen=false;
			dopen=false;
			wp.push_row(2);
			wp(wp.getrnum()-1,0)=i;
			wp(wp.getrnum()-1,1)=j;
			i+=1;j+=1;
		}
		          //when no gap open occur                   //when right gap open,r only needs to compare with s
		else if((ropen==false && dopen==false && r>=d) || (ropen==true && r==start && r>s))          //when down gap open,no right gap is allowed
		{
			ropen=true;
			j+=1;
		}
		else
		{
			dopen=true;
			i+=1;
		}
	}
}
//Menu:[RBLAST]-[_Fillblank]
//Alignment the blank area by using Needleman global alignment technique 
//Take 4 values."_blank" blank area coordiates."_qery" quiry sequence."_bseq" base sequence."alnmatx" standard 0 or 1 waypoint matrix
//Update "alnmatx",which should be formatted into standard 0 or 1 form
void _Fillblank(darray<int> & _blank,sarray<int> & _qeryid,sarray<int> & _bseqid,darray<int> & alnmatx)
{
	int brnum;
	if((brnum=_blank.getrnum())==0)                                             //return if no blank area found
		return;
	else
	{
		darray<int> BKMATX,way;
		sarray<int> subseqa,subseqb;
		for(int b=0;b<brnum;++b)
		{
			if(_blank(b,2)==_blank(b,0) && _blank(b,3)==_blank(b,1))            //if start point equal to end point,no alignment is needed 
			{
				alnmatx(_blank(b,0),_blank(b,1))=1;
				continue;
			}
			else
			{
				subseqa=_qeryid.subseq(_blank(b,0),_blank(b,2));                //intercept sub sequences to prepare for global alignment
				subseqb=_bseqid.subseq(_blank(b,1),_blank(b,3));
				NWalign(subseqa,subseqb,BKMATX);
				_Maxscwp(BKMATX,way);
				for(int i=0;i<way.getrnum();++i)
					alnmatx(_blank(b,0)+way(i,0),_blank(b,1)+way(i,1))=1;       //fill in alnmatx
			}
		}
	}
}

//Menu:[RBLAST]-[_Getscore]
//Compute score for standard 0 or 1 waypoint matrix
//Take 6 values."_scmatx" score matrix."_alnmatx" standard 0 or 1 alignment matrix."_rownum"/"_colnum" row/column number of "_alnmatx"."_qery"/"_bseq" quiry/base sequence."_gappnt" gap penalty
//Return score.
template<class type>
type _Getscore(darray<type> & _scmatx,darray<int> & _alnmatx,const int & _rownum,const int & _colnum,sarray<int> & _qeryid,sarray<int> & _bseqid,const int & _gappnt)
{
	type score=0;
	int gap=0,previ=-1,prevj=-1;
	for(int i=previ+1;i<_rownum;++i)                                           //start from previous point,speed up
	{
		for(int j=prevj+1;j<_colnum;++j)
		{
			if(_alnmatx(i,j)==1)
			{
                if(i!=0 && j!=0)                                               //gaps at the end of alignment were not considered
				   gap=i-previ+j-prevj-2;
				score+=_scmatx(_qeryid[i],_bseqid[j])+(type)gap*_gappnt;     
				previ=i;prevj=j;
				break;
			}
		}
	}
	return score;
}
//Compute score for blank unfilled BLAST alignment matrix
//Take 6 values."_scmatx"scoring matrix."_BLASTmatx" blank unfilled BLAST matrix."_rownum"/"_colnum" row/column number."_qeryid"/"_bseqid" qery and base AA id seq
template<class type>
type _Getscore(darray<type> & _scmatx,darray<int> & _BLASTmatx,const int & _rownum,const int & _colnum,sarray<int> & _qeryid,sarray<int> & _bseqid)
{
	type score=0;int i=0,j=0,prevj=0;
	while(i<_rownum && j<_colnum)                                           
	{
		while(i<_rownum && j<_colnum && _BLASTmatx(i,j)==1)
		{
			score+=_scmatx(_qeryid[i],_bseqid[j]);
			i+=1;j+=1;
			prevj=j;                                                    //always start from the column next to the column where 1 was found.
		}
		j+=1;
		if(j==_colnum)
		{
			i+=1;j=prevj;                                               //i always increases by 1 when one row has been scanned
		}
	}
	return score;
}
		

/*not in use
void AAcomb(const int & _wlen,darray<char> & _AAdb)
{
	int mark1=0,mark2=0,interval;
	_AAdb.fast_resize((int)pow((double)20,_wlen),_wlen);
	for(int w=0; w<_wlen;++w)
	{
		interval=(int)pow((double)20,w);
	    for(int i=0;i<_AAdb.getrnum();++i)
	    {
			if(mark1>19)
				mark1=0;
			_AAdb(i,w)=idaa(mark1);
		    mark2+=1;
			if(mark2==interval)
			{
				mark2=0;
				mark1+=1;
			}
		}
	}
}*/