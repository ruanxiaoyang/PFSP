//////////////////////////////////////////////////////////////
//                         BLOSUM                           //
//            STEVEN HENIKOFF AND JORJA G.HENIKOFF          //
//               Proc. Natl. Acad. Sci. USA                 //
//                Vol. 89, pp. 10915-10919                  //
//           Copyright Belongs to Ruan Xiaoyang             //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
//This module creates blosum scoring matrix from standard block file obtained from http://blocks.fhcrc.org/blocks/uploads/blosum/.
//For each run,a log file was saved at c:\ProFunSite.log indicating the block version and number of blocks contained
//For proper work,each block should be preceded by

//ID   GLU_CARBOXYLATION; BLOCK
//AC   BL00011; distance from previous block=(1,64)
//DE   Vitamin K-dependent carboxylation domain proteins.
//BL   ECA motif; width=40; 99.5%=703; strength=2331

//the sequence should be formatted like 

//ACP2_SPIOL  (    82)  EIGADSLDTVEIVMKLEEEFGVTVEEENAQTITTIQ
//ACP_ECOLI  (     31)  DLGADSLDTVELVMALEEEFDTEIPDEEAEKITTVQ
//ACP_NEUCR  (     86)  DLGLDSLDTVEVVMAIEEEFSIEIPDKDADQIHSVD
//and ended with "//"
//**note:DATA FILE IS "newline sensitive,case sensitive" AND "blank insensitive"

//Protein information was stored in _infvec.block information was stored in _blkfile
//BLOSUM return darray<double> type containing score in bit unit.A detailed report on the obtained blosum matrix was saved at d:\BLOSUM_XX.txt

struct inf_blk                           //structure used to hold information of whole block and titile information
{
	vector<sarray<char>> _checkarr;
	vector<darray<char>> _blkfile;
};

//Menu:[BLOSUM]-[_BLOSUMsingleblock]-[_CprsMatx]
//             -[_Breakstring]
//     [_Impblkfile]-[Breakstring]


//Menu:[_Impblkfile]
//Import block file for BLOSUM use
//Update 1 value."_blkandinf" block and information structure
void _Impblkfile(inf_blk & _blkandinf,const char & _defau)               //import block file from harddisk
{
	clock_t start,end;
	start=clock();
	string fpath;
	char answer,a;
	ifstream ifblkdb;
	bool entseq=false;
	sarray<char> seq;
	darray<char> block;
	string tempstr,ch_a="Blocks";
	vector<sarray<char>> checkarr;
	if(_defau=='Y' ||_defau=='y')
	{
		fpath="c:\\windows\\blocks.dat";
		ifblkdb.open(fpath.c_str());
	}
	else
	{
again:cout<<"Provide path for block database:";
	    cin>>fpath;
		ifblkdb.open(fpath.c_str());
		if(!ifblkdb.is_open())
		{
			cout<<"No file found!\n";
			goto again;
		}
	}
//////////////////////////////////////////////////////////////////////////////////////////check database validity
	ifblkdb.seekg(0,ios::beg);                                                         //set to the start position
    getline(ifblkdb,tempstr);                                                          //obtain the title line of block file
	checkarr=_Breakstring(tempstr);                                                    //format the string line into vector<sarray<char>> type to make the following comparison easier
    if(!(checkarr[0]==ch_a))                                                           //check if the block file has correct title.if not,remind the user
	{
		cout<<"This block file has wrong title,this may(or may not) cause problem.Continue anyway?(Y/N):";
        cin>>answer;	 
        if(answer=='N'||answer=='n')
		{
            ifblkdb.close();
			goto again;
		}
	}
	cout<<"Import block file..."<<endl;
	_blkandinf._checkarr=checkarr;                                                     //assign title information to return value 
	if(_defau=='Y'||_defau=='y')                                                       //to save time,when default block was selected,we do not further import _blkandinf._blkfile
		return;
	while(ifblkdb.get(a))                                                              //search for "ID"--the start marker of a block
	{
		if(a!='I')
			continue;
		else
		{
			ifblkdb.get(a);
			if(a=='D')
				break;
			else
				continue;
		}
	}
	sarray<int> sizeloader;
	int size,n=0;
///////////////////////////////////////////////////////////////////////////////////////import _blkandinf._blkfile
	while(ifblkdb.get(a))
	{
		if(a=='/' || a=='\n')
			continue;
		for(int i=0;i<4;++i)                                                         //skip four lines 
		getline(ifblkdb,tempstr);          
		for(int i=0;i<tempstr.size();++i)
		{
			if(tempstr[i]=='=')
			{
				i+=1;
				while((int)tempstr[i]>=48 && (int)tempstr[i]<=57)
				{
					sizeloader.pushback((int)tempstr[i]-48);
					i+=1;
				}
				break;
			}
		}
		size=0;
		for(int i=0;i<sizeloader.size();++i)
		{
			if(sizeloader[i]!=0)
			   size+=(int)(sizeloader[i]*pow(10.0,sizeloader.size()-i-1));
		}
		sizeloader.clear();
	    entseq=false;
		seq.resize(size);
		while(ifblkdb.get(a))
		{
			if(a=='/')                               
			{
				_blkandinf._blkfile.push_back(block);                                //all blocks were stored in _blkandinf._blkfile
				block.clear();
				break;
			}
			if(entseq==false)
			{
			    if(a==' '|| a=='\n')
				   continue;
		    	if(a==')')  
				   entseq=true;
			}
			else
			{
				if(a!=' ' && a!='\n')
				{
					seq[n]=a;                //store AA in seq until met newline 
				    n+=1;
				}
				else if(a=='\n')
				{
					n=0;
					block.push_row(seq);            //store sequence in block
					entseq=false;
				}
			}
		}
	}
	ifblkdb.close();
	end=clock();
	char time[10],date[10];
	_strtime(time);_strdate(date);
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<" Process time "<<end-start<<" minisec"<<endl;
    outlog<<"Import block file from "<<fpath<<endl;
    for(int i=0;i<checkarr.size();++i)
		outlog<<checkarr[i]<<" ";
	outlog<<"Contain "<<_blkandinf._blkfile.size()<<" Blocks"<<endl<<endl;
}

//menu:[BLOSUM]-[_BLOSUMsingleblock]-[_CprsMatx]
//compress grouping matrix to provide nonredundant grouping information
//take 2 values."_grouparr" upper right nonsymmetrical pairwise correlation matrix of seqs in one block."_rownum" row number of previous matrix
//update "_grouparr" to provide grouping information
void _CprsMatx(darray<int> & _grouparr,const int & _rownum)
{
	for(int i=_rownum-1;i>=1;--i)                                                  //start from lower right corner,search above row for occurance of correlation 
	{
	    for(int m=i-1;m>=0;--m)
		{
	        if(_grouparr(m,i)==1)
			{
				for(int j=i;j<_rownum;++j)
				{
					if(_grouparr(i,j)==1)
					{
						_grouparr(i,j)=0;
						_grouparr(m,j)=1;
					}
				}
				goto nexti;
			}
		}
		for(int j=i+1;j<_rownum;++j)
		{
			if(_grouparr(i,j)==1)
			{
			   for(int m=i-1;m>=0;--m)
			   {
				   if(_grouparr(m,j)==1)
				   {
					  for(int n=i;n<_rownum;++n)
				      {
					      if(_grouparr(i,n)==1)
					      {
						     _grouparr(i,n)=0;
						     _grouparr(m,n)=1;
					     }
				      }
					  goto nexti;
					}
			   }
			}
		}
nexti:;
	}
}

//menu:[BLOSUM]-[_BLOSUMsingleblock]
//calculate mutation frequency count for each block
//take 2 values."_cvgstd" converge standard(also called clustering threshold)."_sngblk" single block.
//update 1 value."Freqcnt" AA mutation frequency count 
int _BLOSUMsingleblock(const int & _cvgstd,darray<char> & _sngblk,darray<double> & Freqcnt)
{
	int rownum=_sngblk.getrnum(),colnum=_sngblk.getcnum(),mark=0,lastrow;
	darray<int> grouparr(rownum,rownum,0),groupinf;
	darray<double> aafreqarr(rownum,3,0.0);
	for(int i=0;i<rownum;++i)                                                         //Calculate percent of identity(sequence similarity)
	{
	   grouparr(i,i)=1;                                                               //for computation convenience,mark every seq similar to itself 
	   for(int k=i+1;k<rownum;++k)
	   {	
	    	mark=0;
		   for(int j=0;j<colnum;++j)
		   {
               if(_sngblk(i,j)==_sngblk(k,j))
                  mark+=1;
		   }
           if((double)mark/(double)colnum*100<(double)_cvgstd)                        //if percent of identity exceed converge standard,mark the two seqs as convergable
			   grouparr(i,k)=0;
		   else
			   grouparr(i,k)=1;
	   }
	}
    _CprsMatx(grouparr,rownum);                                                       //see [_CprsMatx]
	for(int i=0;i<rownum;++i)                                                         //find sequences grouped together,store in darray<int> groupinf
	{
		if(grouparr(i,i)==1) 
		{
			groupinf.push_row(1);
			groupinf(groupinf.getrnum()-1,0)=i;
			for(int j=i+1;j<rownum;++j)
			{
			    	if(grouparr(i,j)>0)
					groupinf.push_to_row(groupinf.getrnum()-1,j);
			}
		}
	}
	if(groupinf.getrnum()<=1)                                                          //if all seqs converged,return 0.which means this block contributed no information
		return 0;
	for(int k=0;k<colnum;++k)                                                          //calculate AA frequency after grouping,store in darray<double> aafreqarr
	{
	   lastrow=0;
	   for(int i=0;i<groupinf.getrnum();++i)
	   {
	    	for(int j=0;j<groupinf.getcnum(i);++j)
		   {                      
			   aafreqarr(lastrow,0)=aaid(_sngblk(groupinf(i,j),k));                    //the 1st column store amino acid ID
			   aafreqarr(lastrow,1)=1.0/groupinf.getcnum(i);                           //the 2nd column store frequency of AA after grouping
			   aafreqarr(lastrow,2)=i;                                                 //the 3rd column store grouping information (for AAs in same group do not contribute to frequency) 
			   lastrow+=1;
		   }
	   }
	   for(int i=0;i<aafreqarr.getrnum()-1;++i)                                        //pairwise mutation count
	   {
		   for(int j=i+1;j<aafreqarr.getrnum();++j)
		   {
			  if(aafreqarr(i,2)==aafreqarr(j,2))                                       //two AAs in same group do not provide information
				   continue;
			  else
			       Freqcnt((int)aafreqarr(i,0),(int)aafreqarr(j,0))+=aafreqarr(i,1)*aafreqarr(j,1);
		   }
	   }
	}
	return 1;
}

//menu:[BLOSUM]
//calculate BLOSUM muation probability matrix from user defined block file
//take 3 values."_cvgstd" converge standard."_logbase" base number of log."_blkandinf" a structure containing block file and associated version number
//return "_logbase" based score matrix
darray<double> BLOSUM(const int & _cvgstd,const int & _logbase,const char & _defau='Y')
{
	inf_blk _blkandinf;
	_Impblkfile(_blkandinf,_defau);
	bool found=false;
	char a,answer;int b;string c; double RF;                                                           //"RF" short for relative frequency
	vector<sarray<char>> checkarr;                                                              //hold block version information
    darray<double> Sij_log(20,20,0.0);
	ifstream read("c:\\windows\\_ProFunSit_BLOSUM_RltvFreq");                                   //look up existing record to see if logSij can be directly calculated from known RltvFreq
    if(read)
	{
		while(read.get(a))
		{
			if(a!='>')                                                                          //search for start marker '>'
				continue;
			read>>b;
			if(b!=_cvgstd)                                                                      //if converge standard do not match,this is not the one we need
				continue;
			getline(read,c);
			checkarr=_Breakstring(c);                                                           //format the string line into vector<sarray<char>> type to make the following comparison easier
			if(checkarr.size()!=_blkandinf._checkarr.size())
				continue;
			for(int i=0;i<checkarr.size();++i)
			{
				if((checkarr[i]==_blkandinf._checkarr[i])==false)                               //give up whenever matching failed
					goto nexta;
			}
			for(int i=0;i<20;++i)                                                               //if the imported block file have same label with known record,directly use record data and return
			{
			    for(int j=0;j<20;++j)
			    {
			        read>>RF;
					if(j>=i)
			        Sij_log(j,i)=Sij_log(i,j)=10*log(RF)/log((double)_logbase);                           //calculate "_logbase" based logarithm of relative frequency
			    }
			}
			if(_defau=='N' || _defau=='n')
			{
				found=true;
				cout<<"Existing record found with same block file name,converge standard and log base\n"
					<<"Still output result?[Y/N]:";
				cin>>answer;
				if(answer=='Y'||answer=='y')
					break;
			}
			return Sij_log;
nexta:;	
		}
	}
	////////////////////////////////////////////////////////////////////////////////////////////
	if(_defau=='Y'||_defau=='y')                                                              //if use default block file,only _blkandinf._checkarr was filled by _Impblkfile.(_blkandinf._blkfile was not filled in order to save time)
	{                                                                                         //this can prevent 1:default was used but no corresponding record was found in c:\\windows\\_ProFunSit_BLOSUM_RltvFreq
		_blkandinf._checkarr.clear();                                                         //                 2:when non-default was used,we need not to import again
		_Impblkfile(_blkandinf,'N');
	}
    darray<double> Freqcnt(20,20,0.0),Freq(20,20,0.0),AA_AAProb(20,20,0.0),RltvFreq(20,20,0.0); //"Freqcnt" frequency count."Freq"="Freqcnt"/total frequency count."AA_AAProb" random frequency
	sarray<double> AAFreq(20,0.0);
	double TotalFreqcnt=0,temp=0;
	int blocknum=0,seqnum=0,colnum=0,AAnum=0;
	cout<<"Constructing block...";
    for(int i=0;i<_blkandinf._blkfile.size();++i)
	{
		PctMarker(i,_blkandinf._blkfile.size()-1,2);
		blocknum+=_BLOSUMsingleblock(_cvgstd,_blkandinf._blkfile[i],Freqcnt);
		seqnum+=_blkandinf._blkfile[i].getrnum();                                               //record number of sequence in each block
		colnum+=_blkandinf._blkfile[i].getcnum(0);                                              //record number of column in each block
		AAnum+=_blkandinf._blkfile[i].getrnum()*_blkandinf._blkfile[i].getcnum(0);              //record total AA number in each block
	}
	cout<<endl;
    for(int i=0;i<20;++i)
	{
		for(int j=i;j<20;++j)
		{
			if(j!=i)
			{
				Freqcnt(i,j)+=Freqcnt(j,i);                        //add A->B to B->A mutation 
			    Freqcnt(j,i)=Freqcnt(i,j);    
				if(Freqcnt(i,j)==0)
				{
					cout<<"0 count occurred,re-clustering fail!\n";
					return Sij_log;
				}
				TotalFreqcnt+=Freqcnt(i,j);                        //add up total frequency count
			}
			if(j==i)
				TotalFreqcnt+=Freqcnt(i,i);
		}
	}
    for(int i=0;i<20;++i)
	{
		for(int j=0;j<20;++j)
		{
            if(i<=j)
			   Freq(i,j)=Freqcnt(i,j)/TotalFreqcnt;
			if(i==j)
                temp+=Freqcnt(i,j);                                //temp was used to calculate AA frequency
			else
				temp+=Freqcnt(i,j)/2.0;
		}
		AAFreq[i]=temp/TotalFreqcnt;                               //calculate AA frequency
		temp=0;
	}
	for(int i=0;i<20;++i)
	{
		for(int j=i;j<20;++j)
		{
			if(i==j)
   			    AA_AAProb(i,j)=AAFreq[i]*AAFreq[j];
			else
				AA_AAProb(i,j)=2*AAFreq[i]*AAFreq[j];
			RltvFreq(i,j)=Freq(i,j)/AA_AAProb(i,j);
			Sij_log(i,j)=10*log(RltvFreq(i,j))/log((double)_logbase);
			if(i!=j)
				Sij_log(j,i)=Sij_log(i,j);                          //make "Sij_log" symmetric
		}
	}
	if(found==false)
	{
		ofstream recRF("c:\\windows\\_ProFunSit_BLOSUM_RltvFreq",ios_base::app);                    //record converge standard,block file information and relative frequency matrix
		recRF<<endl<<">"<<_cvgstd<<" ";
		for(int i=0;i<_blkandinf._checkarr.size();++i)
			recRF<<_blkandinf._checkarr[i]<<" ";
		recRF<<endl;
		RltvFreq.record(recRF);
	}
//////////////////////////////////////////////////////////////////////////////////////////////////////output detailed result
	MKDIR("C:\\PFSP\\BLOSUM");
	string path;
	path="C:\\PFSP\\BLOSUM\\BLOSUM_"+inttostr(_cvgstd)+".txt";
	ofstream outBLS(path.c_str());
	char time[10],date[10];_strtime(time);_strdate(date);
	outBLS<<date<<" "<<time<<endl;
	outBLS<<"Protein Functional Sites Prediction--BLOSUM:(C++) Copyright 2009,Ruan Xiaoyang"<<endl;
	for(int i=0;i<_blkandinf._checkarr.size();++i)
		outBLS<<_blkandinf._checkarr[i]<<" ";
    outBLS<<endl<<"Re-clustering percentage "<<_cvgstd<<"%"<<endl<<endl;
	outBLS<<_blkandinf._blkfile.size()<<" blocks processed  "<<blocknum<<" blocks contributed pairs to matrix"<<endl;
	outBLS<<TotalFreqcnt<<" total pairs  "<<seqnum<<" total sequences  "<<colnum<<" total columns  "<<AAnum<<" total AAs"<<endl;
	outBLS<<"AA frequencies"<<endl;
	for(int i=0;i<AAFreq.size();++i)
		outBLS<<idaa(i)<<"         ";
	outBLS<<endl;
	outBLS<<AAFreq<<endl<<endl;

	outBLS<<"Number of pairs count"<<endl;
    for(int i=0;i<20;++i)
	    outBLS<<idaa(i)<<"         ";
	outBLS<<endl; 
	outBLS<<Freqcnt;

    outBLS<<"Frequency"<<endl;
    for(int i=0;i<20;++i)
	    outBLS<<idaa(i)<<"         ";
	outBLS<<endl;
	outBLS<<Freq;

	outBLS<<"Relative odds"<<endl;
	for(int i=0;i<20;++i)
	    outBLS<<idaa(i)<<"         ";
	outBLS<<endl;
	outBLS<<RltvFreq;

	outBLS<<"Sij_log "<<_logbase<<endl;
	for(int i=0;i<AAFreq.size();++i)
	    outBLS<<idaa(i)<<"         ";
	outBLS<<endl;
	outBLS<<Sij_log;
	outBLS.close();
	TMPPATH=path;
	DSPFILE(TMPPATH);
	return Sij_log;
}




