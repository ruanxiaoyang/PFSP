//////////////////////////////////////////////////////////////
//             Multiple Sequences Alignment                 //
//           Copyright Belongs to Ruan Xiaoyang             //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
//This module perform fast and accurate MSA.Several modes of alignment is available.The work flow including genetic distance estimation (accurate or approximate),evolutionary tree construction and multiple sequence alignment.
//For each run,a log file was saved at c:\ProFunSite.log indicating the alignment information

//menu: [MSA]-[BysDistMatx]-[BayesianDist]
//           -[AprxDistMatx]-[AprxDist]
//           -[NJtree]
//           -[_Alignseqs]-[_SeqMge]-[Length]
//                        -[_Iterate]
//                        -[_OmitBlank]
//           -[_ImpDistMatx]
//           -[_MSAsc]
//           -[_Iterate]
//           -[_MdfMSA]

//Menu:[MSA]-[_ImpDistMatx]
//Import distance and variance matrix file from record
//Take 1 value."_imp" ifstream
//Output 2 values."distmatx" distance matrix."varmatx" variance matrix
void _ImpDistMatx(ifstream & _imp,darray<double> & distmatx,darray<double> & varmatx)
{
	int seqnum;
	double var;
	_imp>>seqnum;
	distmatx.fast_resize(seqnum,seqnum);
	varmatx.fast_resize(seqnum,seqnum);
	for(int i=0;i<seqnum;++i)
	{
		for(int j=0;j<seqnum;++j)
		{
			_imp>>var;
			distmatx(i,j)=var;
		}
	}
	for(int i=0;i<seqnum;++i)
	{
		for(int j=0;j<seqnum;++j)
		{
			_imp>>var;
			varmatx(i,j)=var;
		}
	}
	_imp.close();
}
//Menu:[MSA]-[_MSAsc]
//Calculate score for MSA.The scoring matrix is Gonnet
//Take 2 value."_seqsidmatx"."_gappnt" gap penalty
//Return score;
double _MSAsc(darray<int> & _seqsidmatx,const int & _gappnt=-3)
{
	int I,J,row=_seqsidmatx.getrnum();
	double sc=0.0;
	for(int j=0;j<_seqsidmatx.getcnum();++j)
	{
		for(int i=0;i<row;++i)
		{
			for(int ii=i+1;ii<row;++ii)
			{
				I=_seqsidmatx(i,j);
				J=_seqsidmatx(ii,j);
				if(I>=0 && J>=0)
				   sc+=Gonnet(I,J);
				else if ((I<0 && J>=0)||(I>=0 && J<0))
					sc+=(double)_gappnt;
			}
		}
	}
	sc/=((double)(row*(row-1))/2.0);
	return sc;
}

//Menu:[MSA]-[_Iterate] [MSA]-[_Alignseqs]-[_Iterate]
//This module refines the MSA result to correct obvious errors
//Take 3 values."_seqsidmatx" MSA resulted aligned AA id sequences."_ancthld" anchoring point threshold."_amp"gap amplification factor
//Return 1 value."posicnt"position AA count of the rectified matrix.Update 1 value "_seqsidmatx"
void _Iterate(darray<int> & _seqsidmatx,const int & _ancthld,const int & _amp,darray<int> & posicnt)
{
	int row=_seqsidmatx.getrnum(),col=_seqsidmatx.getcnum(),_index_,J,origsc,origdist,markJ,aaid;
	posicnt.fast_resize(20,col,0);
	sarray<int> anchor,ancdist(col);
	for(int j=0;j<col;++j)
	{
		for(int i=0;i<row;++i)
		{
			if(_seqsidmatx(i,j)>=0)
				posicnt(_seqsidmatx(i,j),j)+=1;
		}
		if(posicnt.columnmax(j,0,_index_)>=(int)((double)(row*_ancthld)/100.0))
			anchor.pushback(j);
	}
	if(anchor.size()==0)
	{
        anchor.resize(col);
		for(int j=0;j<col;++j)
			anchor[j]=j;
	}
	int a=0,pnt=_amp*((int)(row/10.0)==0?1:(int)(row/10.0));                                           //calculate the min distance of each position to anchoring point.
    for(int j=0;j<col;++j)                                                                             //multiply the distance with amplification factor
	{
		if(a!=anchor.size() && j<anchor[a])
		{
			if(a==0)
			   ancdist[j]=pnt*(anchor[a]-j);
			else if(a>0 && a<=anchor.size()-1)
				ancdist[j]=((j-anchor[a-1])>(anchor[a]-j)?pnt*(anchor[a]-j):pnt*(j-anchor[a-1]));
		}
		else if(a!=anchor.size() && j==anchor[a])
		{
			ancdist[j]=0;
			a+=1;
		}
		else if(a==anchor.size())
			ancdist[j]=pnt*(j-anchor[a-1]);
	}
	bool shifted=false;                                                         //"shifted" controls whether a new round of iteration is needed
	for(int i=0;i<row;++i)
	{
		for(int j=0;j<col;++j)
		{
			if(_seqsidmatx(i,j)<0)
				continue;
			else
			{
				aaid=_seqsidmatx(i,j);
				_seqsidmatx(i,j)=-4;
			    origdist=ancdist[j];
				origsc=posicnt(aaid,j)-origdist;                                //score of a position is the count number of the AA minus gap penalty  
				J=j;
				while(J>=0 && _seqsidmatx(i,J)<0)                               //shift to the left gap
					J-=1;
				J+=1;
				while(J<col && _seqsidmatx(i,J)<0)                              //search the blank space for a better position
				{
					if((posicnt(aaid,J)-ancdist[J])>origsc || ((posicnt(aaid,J)-ancdist[J])==origsc && ancdist[J]<=origdist))
					{     //J has a better score                                  //J has same score but lower distance                               
						markJ=J;
						origdist=ancdist[J];
						origsc=posicnt(aaid,J)-origdist;
					}
					J+=1;
				}
				_seqsidmatx(i,markJ)=aaid;
				if(markJ!=j)
				{
					shifted=true;
					posicnt(aaid,markJ)+=1;                                      //update the position count
					posicnt(aaid,j)-=1;
				}
				j=J-1;
			}
		}
	}
	if(shifted==true)                                                            //if shift was detected,run another round of iteration untill convergence
		_Iterate(_seqsidmatx,_ancthld,_amp,posicnt);
	return;
}

//Menu:[MSA]-[_OmitBlank] [MSA]-[_Alignseqs]-[_OmitBlank]
//Omit postions that was purely composed of gaps
//Take 2 values."_seqsICmatx" int or char seqs."_posicnt" position count
//Output 1 value."identical" sarray storing positions with identical AA
//Update 1 value."_seqsICmatx" pure gap will be omitted
template <class type>
void _OmitBlank(darray<type> & _seqsICmatx,darray<int> & _posicnt,sarray<char> & identical)
{
	int colnum=0,max,_index_;
	sarray<int> nonemptyposi;
	darray<type> tmpseqsICmatx;
	identical.clear();
	for(int j=0;j<_posicnt.getcnum();++j)
	{
		if((max=_posicnt.columnmax(j,0,_index_))!=0)
		{
			colnum+=1;
			nonemptyposi.pushback(j);
			identical.pushback(max==_seqsICmatx.getrnum()?'*':' ');
		}
	}
	if(colnum==_seqsICmatx.getcnum())                                          //if no pure gap found,no update is needed
	   return;
	else
	{
		tmpseqsICmatx.fast_resize(_seqsICmatx.getrnum(),colnum);
		for(int j=0;j<colnum;++j)
		{
			for(int i=0;i<_seqsICmatx.getrnum();++i)
				tmpseqsICmatx(i,j)=_seqsICmatx(i,nonemptyposi[j]);
		}
		_seqsICmatx=tmpseqsICmatx;
	}
}
//Menu:[MSA]-[_MdfMSA]
//Modify final MSA result.This module will change "_seqsidmatx"
//Take 1 value."_seqsidmatx"AA id sequences matrix
//Output 4 values."dspmatx"display matrix."identical"sarray storing positions with identical AA."ancthld"anchoring point threshold."amp"gap amplification factor
double _MdfMSA(darray<int> & _seqsidmatx,darray<char> & dspmatx,sarray<char> & identical,int & ancthld,int & amp,const int & _gappnt=-3)
{
	char answer;
    int trial=0;                                                                                 //control output file name
	double score;
	darray<int> tmpmatx,posicnt;
again:	cout<<"Specify [anchoring point threshold][gap amplification factor] seperate by space\n"
		<<"ep:90 3 represent [anchoring point threshold 90%][gap amplification factor 3]\n"
		<<"Input:";
	cin>>ancthld>>amp;
	trial+=1;
	tmpmatx=_seqsidmatx;                                                                         //keep original "_seqsidmatx" intact
	cout<<"Iterating...\n";
	_Iterate(tmpmatx,ancthld,amp,posicnt);                                                       //iterate final MSA result
    _OmitBlank(tmpmatx,posicnt,identical);                                                       //omit pure gap area(blank) after iteration
	score=_MSAsc(tmpmatx,_gappnt);                                                                       //calculate score of MSA after iteration
	MatxItoC(tmpmatx,dspmatx);                                                                   //convert AA id seqs to AA char seqs
	string prvwpath="c:\\PFSP\\TEMP\\Final_Trial_"+inttostr(trial)+"_MSA_Preview_"+inttostr(ancthld)+"_"+inttostr(amp)+".txt"; //preview file path
	ofstream prvwout(prvwpath.c_str());                                                          //preview file ofstream
	prvwout<<"This MSA result is for your preview only.Please delete this temp file in C:\\PFSP\\TEMP\n"
		<<"Tips:Increase gap amplification factor to reduce gap.Use moderately lower anchor threshold to allow for more anchoring point\n" 
		<<"Anchor threshold:"<<ancthld<<"  Gap amplification factor:"<<amp<<"  Length:"<<dspmatx.getcnum()<<"  Length diff:"<<_seqsidmatx.getcnum()-dspmatx.getcnum()
		<<"  Score(Gonnet):"<<score<<endl;
	LandMarker(dspmatx.getcnum(),prvwout);prvwout<<endl;
	MileStone(dspmatx.getcnum(),prvwout);prvwout<<endl;
	prvwout<<identical<<endl;
	prvwout<<dspmatx;
	prvwout<<identical<<endl;
	MileStone(dspmatx.getcnum(),prvwout);prvwout<<endl;
	LandMarker(dspmatx.getcnum(),prvwout);prvwout<<endl;
	prvwout.close();
	TMPPATH=prvwpath;                                                                             
	DSPFILE(TMPPATH);
	cout<<"Are you satisfied with this MSA result?[Y/N]:";
	cin>>answer;
	if(answer=='N' || answer=='n')
		goto again;
	else
        _seqsidmatx=tmpmatx;
	return score;
}

//Menu:[_MdfMSA]
//Overloaded function of [MSA]-[_MdfMSA].Specifically designed for mainmenu option i.It refines previous MSA record
//Take 1 value."_dspmatx" AA char matrix(usually imported from MSA record file)
void _MdfMSA(darray<char> & _dspmatx)
{
	int ancthld,amp;
	darray<char> dspmatx;
	darray<int> seqsidmatx;
	sarray<char> identical;
	MatxCtoI(_dspmatx,seqsidmatx);
	_MdfMSA(seqsidmatx,dspmatx,identical,ancthld,amp);
	MKDIR("C:\\PFSP\\MSA");
	string record;
	record="C:\\PFSP\\MSA\\MSA_"+inttostr(dspmatx.getrnum())+"_Refined_"+inttostr(ancthld)+"_"+inttostr(amp)+"_record";
    ofstream os(record.c_str());
	os<<dspmatx.getrnum()<<" "<<dspmatx.getcnum()<<endl;
	os<<identical<<endl;
	os<<dspmatx;
	cout<<"Refined MSA has been saved to "<<record<<endl;
}

//Menu:[MSA]-[_Alignseqs]
//This module align multiple seqs.It is the main component of [MSA] 
//"_nbtree" structure hold evolutionary tree(see EvoTree for detail)
//"_mwit"midway iteration switcher."_mancthld"midway anchoring point threshold."_mamp"midway gap amplification factor
//Other parameters see [MSA]
template<class type>
void _Alignseqs(vector<sarray<int>> & _multiarrid,nbstr & _nbtree,darray<double> & _distmatx,vector<darray<type>> & _scoredb,darray<int> & seqsmatx,sarray<char> & identical,const int & _wlen=3,const int & _gappnt=-3,const int & _cutoffper=95,const char & _DorB='D',const char & _mwit='N',const int & _mancthld=95,const int & _mamp=3)
{
	cout<<"Align multiple sequences...";
	int dist,qeryid,bseqid;
	type tempscore;
	vector<darray<int>> alignarrvec;
	darray<int> combseq,alignmatx,posicnt;
	for(int i=0;i<_nbtree._nbsize;++i)
	{
		PctMarker(i,_nbtree._nbsize-1);
		if(_nbtree._lkbk(i,0)<0 && _nbtree._lkbk(i,1)<0)                                                        //when two seqs are both single seq,use RBLAST alignment
		{
			qeryid=_nbtree._marr1(i,0);
			bseqid=_nbtree._marr2(i,0);
			int tmp=(int)(_distmatx(qeryid,bseqid)*100);
			dist=(tmp<1?1:tmp);
			RBLAST(_wlen,_scoredb[dist-1],_multiarrid[qeryid],_multiarrid[bseqid],dist,alignmatx,tempscore,_gappnt,_cutoffper,_DorB);
			_SeqMge(_multiarrid[qeryid],_multiarrid[bseqid],alignmatx,combseq);                                 //create aligned seqs
		}
		else if(_nbtree._lkbk(i,0)>=0 && _nbtree._lkbk(i,1)<0)                                                  //when one part is aligned seqs,use global alignment rule
			_SeqMge(alignarrvec[_nbtree._lkbk(i,0)],_multiarrid[_nbtree._marr2(i,0)],combseq,_gappnt);
		else if(_nbtree._lkbk(i,0)<0 && _nbtree._lkbk(i,1)>=0)
			_SeqMge(_multiarrid[_nbtree._marr1(i,0)],alignarrvec[_nbtree._lkbk(i,1)],combseq,_gappnt);
		else if(_nbtree._lkbk(i,0)>=0 && _nbtree._lkbk(i,1)>=0)
			_SeqMge(alignarrvec[_nbtree._lkbk(i,0)],alignarrvec[_nbtree._lkbk(i,1)],combseq,_gappnt);
		if(combseq.getrnum()>=3 && (_mwit=='Y'||_mwit=='y'))                                                    //do midway iteration if it is enabled,
		{
			_Iterate(combseq,_mancthld,_mamp,posicnt);
            _OmitBlank(combseq,posicnt,identical);
		}
		alignarrvec.push_back(combseq);
		combseq.clear();
	}
	seqsmatx=alignarrvec[alignarrvec.size()-1];
	if(identical.size()==0)
	{
		int I;
		identical.resize(seqsmatx.getcnum());
		for(int j=0;j<seqsmatx.getcnum();++j)
		{
			for(I=1;I<seqsmatx.getrnum();++I)
			{
				if(seqsmatx(I,j)!=seqsmatx(0,j))
					break;
			}
			identical[j]=(I==seqsmatx.getrnum()?'*':' ');
		}
	}
	alignarrvec.clear();
	cout<<endl;
}

//menu:[MSA]-[_Alignseqs]-[_SeqMge]
//merge two sarray<char> sequence according to "_matrix" to give two properly aligned seqs stored in darray<char>
//take 3 values."_qery"/"_bseq" query and base sequence."_matrix" standard 0 or 1 alignment matrix(row for "_qery",column for "_bseq").
//output 1 value."combseq" two aligned seqs stored in darray<char>,first row store "_bseq",second row store "_qery",gap was represented by '-' char
void _SeqMge(sarray<int> & _qeryid,sarray<int> & _bseqid,darray<int> & _matrix,darray<int> & combseq)
{
	int rownum=_qeryid.size(),colnum=_bseqid.size(),previ=-1,prevj=-1,posi=-1,size,remain,gapi,gapj;
    size=Length(_matrix,rownum,colnum);                                 //calculate total length of aligned seq,use it to initiate "combseq"
	combseq.fast_resize(2,size);                                         //specify volume of "combseq"
	for(int i=0;i<rownum;++i)
	{
		for(int j=prevj+1;j<colnum;++j)
		{
			if(_matrix(i,j)==1)
			{
				gapi=i-previ-1;                                          //count gap
				gapj=j-prevj-1;
				if(gapi>0)
				{
					for(int g=1;g<=gapi;++g)
					{
						combseq(0,posi+g)=-4;                          //"posi" recorded current position in "combseq" 
                        combseq(1,posi+g)=_qeryid[previ+g];
					}
					posi+=gapi;
				}
				if(gapj>0)
				{
					for(int g=1;g<=gapj;++g)
					{
						combseq(0,posi+g)=_bseqid[prevj+g];
						combseq(1,posi+g)=-4;
					}
					posi+=gapj;
				}
				posi+=1;
				previ=i;prevj=j;
				combseq(0,posi)=_bseqid[j];
				combseq(1,posi)=_qeryid[i];
				if(prevj==colnum-1)                                      //if reach the end of colunm,check if row gap remain,update "combseq" and return
				{
					if((remain=rownum-1-i)>0)
					{
						for(int g=1;g<=remain;++g)
						{
							combseq(0,posi+g)=-4;
							combseq(1,posi+g)=_qeryid[previ+g];
						}
					}
					return;
				}
				break;
			}
		}
	}
	if((remain=colnum-1-prevj)>0)                                        //if reach the end of row,check if column gap remain,update "combseq" and return
	{
	    for(int g=1;g<=remain;++g)
		{
			combseq(0,posi+g)=_bseqid[prevj+g];
			combseq(1,posi+g)=-4;
		} 
	}

	return;
}
//Menu:[MSA]-[_Alignseqs]-[_SeqMge]
//Merge two darray<char> seqs into a combined one according to Needleman global rule
//Take 3 values."_marr1"/"_marr2" multiple(or sinlge) seqs store as darray<char>."_gappnt" gap penalty
//Output 1 value."combseq" final combination of former seqs 
void _SeqMge(darray<int> & _marrid1,darray<int> & _marrid2,darray<int> & combseq,const int & _gappnt)
{
	int rownum=_marrid1.getcnum(0),colnum=_marrid2.getcnum(0),sn1=_marrid1.getrnum(),sn2=_marrid2.getrnum(),id1,id2,gapi,gapj,previ=-1,prevj=-1,posi=-1,size,remain;
    double tmpscore,s,max;
	darray<int> stdmatx(rownum,colnum,0);
	darray<double> BKMATX(rownum+1,colnum+1,0.0);
	for(int j=colnum-1;j>=0;--j)                                                   //calculate BKMATX
	{
		for(int i=rownum-1;i>=0;--i)
	    {		
		    tmpscore=0;
			for(int k=0;k<sn2;++k)
			{
				for(int m=0;m<sn1;++m)
				{
					id1=_marrid1(m,i);
					id2=_marrid2(k,j);
					if(id1>=0 && id2>=0)
				       tmpscore+=Gonnet(id1,id2);
					else if(id1<0 && id2<0)
						tmpscore+=2;
				}
			}
			s=tmpscore+BKMATX(i+1,j+1);   
			max=(BKMATX(i+1,j)>BKMATX(i,j+1)?BKMATX(i+1,j):BKMATX(i,j+1));
			BKMATX(i,j)=(s>max?s:max);
		} 
	}
	int I=0,J=0,rgap=0,dgap=0;
	double rsc,dsc,ssc,r,d,start;
    while(I<rownum && J<colnum)                                                     //calculate standard 0 or 1 waypoint matrix
	{
		start=BKMATX(I,J);
		rgap+=1;
	    dgap+=1;
		r=BKMATX(I,J+1);
		d=BKMATX(I+1,J);
		s=BKMATX(I+1,J+1);
		rsc=r+rgap*_gappnt;
		dsc=d+dgap*_gappnt;
		ssc=s;
        if(((r!=start && d!=start) || (s>=r && s>=d)||(ssc>=rsc && ssc>=dsc))|| (rgap>1 && (r!=start || s>=r || ssc>=rsc)) || (dgap>1 && (d!=start || s>=d|| ssc>=dsc)))
		{
			rgap=0;dgap=0;
			stdmatx(I,J)=1;
			I+=1;J+=1;
		}
		else if((rgap==1 && dgap==1 && (r>d ||rsc>=dsc)) || (rgap>1 && r>s && rsc>ssc))
		{
			dgap=0;
			J+=1;
		}
		else
		{
			rgap=0;
			I+=1;
		}
	}
    size=Length(stdmatx,rownum,colnum);                                       //calculate the volume of combined sequence
	combseq.fast_resize(sn1+sn2,size);                                        
    for(int i=0;i<rownum;++i)                                                 //merge the two multiarray into a single one according to above alignemnt matrix
	{
		for(int j=prevj+1;j<colnum;++j)
		{
			if(stdmatx(i,j)==1)
			{
				gapi=i-previ-1;
				gapj=j-prevj-1;
				if(gapi>0)
				{
					for(int g=1;g<=gapi;++g)
					{
						for(int sn=0;sn<sn2;++sn)
						    combseq(sn,posi+g)=-4;
						for(int sn=0;sn<sn1;++sn)
                            combseq(sn2+sn,posi+g)=_marrid1(sn,previ+g);
					}
					posi+=gapi;
				}
				if(gapj>0)
				{
					for(int g=1;g<=gapj;++g)
					{
						for(int sn=0;sn<sn2;++sn)
						    combseq(sn,posi+g)=_marrid2(sn,prevj+g);
                        for(int sn=0;sn<sn1;++sn)
						    combseq(sn2+sn,posi+g)=-4;
					}
					posi+=gapj;
				}
				posi+=1;
				previ=i;prevj=j;
                for(int sn=0;sn<sn2;++sn)
				    combseq(sn,posi)=_marrid2(sn,j);
                for(int sn=0;sn<sn1;++sn)
				    combseq(sn2+sn,posi)=_marrid1(sn,i);
				if(prevj==colnum-1)
				{
					if((remain=rownum-1-i)>0)
					{
						for(int g=1;g<=remain;++g)
						{
                            for(int sn=0;sn<sn2;++sn)
							    combseq(sn,posi+g)=-4;
                            for(int sn=0;sn<sn1;++sn)
							    combseq(sn2+sn,posi+g)=_marrid1(sn,previ+g);
						}
					}
					return;
				}
				break;
			}
		}
	}
	if((remain=colnum-1-prevj)>0)
	{
	    for(int g=1;g<=remain;++g)
		{
			for(int sn=0;sn<sn2;++sn)
			    combseq(sn,posi+g)=_marrid2(sn,prevj+g);
            for(int sn=0;sn<sn1;++sn)
			    combseq(sn2+sn,posi+g)=-4;
		} 
	}
	return;
}
//overloaded functon of _SeqMge(darray<char> & _marr1,darray<char> & _marr2,darray<char> & combseq) to take sarray<char> type in the first argument
void _SeqMge(sarray<int> & _arrid1,darray<int> & _marrid2,darray<int> & combseq,const int & _gappnt)
{
	darray<int> marrid1;
	marrid1.push_row(_arrid1);
	_SeqMge(marrid1,_marrid2,combseq,_gappnt);
}
//overloaded functon of _SeqMge(darray<char> & _marr1,darray<char> & _marr2,darray<char> & combseq) to take sarray<char> type in the second argument
void _SeqMge(darray<int> & _marrid1,sarray<int> & _arrid2,darray<int> & combseq,const int & _gappnt)
{
	darray<int> marrid2;
	marrid2.push_row(_arrid2);
	_SeqMge(_marrid1,marrid2,combseq,_gappnt);
}


//This module takes 11 values.
//"_seqandinf" structure hold sequence file and its information.such as "fastaformat" 
//"_scoredb" Dayhoff PAM scoring matrices
//"_EXorAP" choose 'E/e'or 'A/a' for bayesian exact distance or approximate distance
//"_distype" when "_EXorAP" was 'E/e',choose 'B/b' for blast local distance.'N/n'for needleman global distance.when"_EXorAP" was 'A/a',choose 'S/s'for simple p distance,'P/p'for poisson distance,'G/g'for gamma distance
//"_wlen" word length of library (see RBLAST for detail)
//"_gappnt" gap penalty
//"_cutoffper" cutoff percent.the percentile threshold up which a hit and a maximum alignment will be considered significant(see RBLAST for detail)
//"_DorB" amino acid frequency data used to generate random AA sequence.'D'for dayhoff,'B'for blosum(see PAM for detail)
//"_gma" gamma parameter when use gamma approximate distance(see EvoDist for detail)
//"_step" control bayesian accuracy level(see BayesianDist for detail)
//"_acculvl" control bayesian accuracy level(see BayesianDist for detail)
//Output 2 values
//"seqsidmatx" aligned sequence id matrix
//"dspmatx" aligned sequence matrix
//"identical" a char array holds identical AA column infor
void MSA(proformat & _seqsandinf,vector<darray<double>> & _scoredb,darray<int> & seqsidmatx,darray<char> & dspmatx,sarray<char> & identical,const char & _EXorAP,const char & _distype,const int & _wlen=3,const int & _gappnt=-3,const int & _ancthld=50,const int & _amp=2,const int & _cutoffper=95,const char & _DorB='D',const int & _gma=1,const int & _step=10,const int &_acculvl=4)
{
	char answer;
	double score=(double)INT_MIN;
	darray<double> distmatx,varmatx;                                                                             //distance matrix,variance matrix
	darray<int> seqsmatxuns;                                                                                     //unsorted AAid seqs alignment matrix
	sarray<int> seqsidorder;                                                                                     //sequences id order,used to restore seq order disrupted by MSA
	string distmatxrecpath;
	answer='N';
	MKDIR("c:\\PFSP\\DistMatx");                                                                                 //create file folder to hold distance matrix record
	if(_EXorAP =='E' || _EXorAP=='e')                                                                       
       distmatxrecpath="c:\\PFSP\\DistMatx\\DistMatx_"+inttostr(_seqsandinf._seqnum)+"_"+string(&_EXorAP,sizeof(char))+"_"+string(&_distype,sizeof(char))+"_"+inttostr(_wlen)+"_"+inttostr(_step)+"_"+inttostr(_acculvl);
	else if(_EXorAP =='A' || _EXorAP=='a')
	   distmatxrecpath="c:\\PFSP\\DistMatx\\DistMatx_"+inttostr(_seqsandinf._seqnum)+"_"+string(&_EXorAP,sizeof(char))+"_"+string(&_distype,sizeof(char))+"_"+inttostr(_gma);
	ifstream impdm(distmatxrecpath.c_str());
	if(impdm)
	{
		cout<<"Same distance matrix record found!Use this record file?[Y/N]:";                                   //check existing record
		cin>>answer;
		if(answer=='Y'||answer=='y')
			_ImpDistMatx(impdm,distmatx,varmatx);                                                                //import record distance and variance matrix if user agreed."impdm" was closed in "_ImpDistMatx"function.
	}
	if(!impdm || answer=='N'||answer=='n')                                                                       //if no record found or user do not agree to use previous record
	{
		if(_EXorAP =='E' || _EXorAP=='e')
		{
			cout<<"Creating pair wise bayesian exact distance matrix...";
			BysDistMatx(_seqsandinf._multiseqid,distmatx,varmatx,_scoredb,_distype,_wlen,_gappnt,_step,_acculvl);//obtain beiyesian distance matrix
		}
		else if(_EXorAP =='A' || _EXorAP=='a')
		{
			cout<<"Creating pair wise approximate distance matrix...";
			AprxDistMatx(_seqsandinf._multiseqid,distmatx,varmatx,_distype,_gma);                                //obtain approximate distance matrix
		}
		ofstream outdm(distmatxrecpath.c_str());                                                                 //write distance and variance record to specific directory
		outdm<<_seqsandinf._seqnum<<endl;
		outdm<<distmatx;
		outdm<<varmatx;
		outdm.close();
	}
	nbstr nbmatx;                                                                                                //structure hold genetic tree information(see "EvoTree" for detail)
	cout<<endl<<"Constructing NJ tree...\n";
	NJtree(distmatx,nbmatx);                                                                                     //create neighbor joining unrooted tree
	seqsidorder=nbmatx._marr1[nbmatx._nbsize-1]+nbmatx._marr2[nbmatx._nbsize-1];                                 //obtain seqs id order of the aligned seqs(this is used to restore the original order appeared in user's file)
	
	int mancthld=INT_MIN,mamp=INT_MIN,mtrial=0;                                                                  //midway anchoring point threshold,midway gap amplification factor
	cout<<endl<<"Enable midway iteration?[Y/N]:";                                                                //midway iteration(also called shifting,the multiple alignment result will be iterated each time more than two sequences were aligned).
	cin>>answer;
	if(answer=='N'||answer=='n')
		_Alignseqs(_seqsandinf._multiseqid,nbmatx,distmatx,_scoredb,seqsmatxuns,identical,_wlen,_gappnt,_cutoffper,_DorB);
	else 
	{
		MKDIR("c:\\PFSP\\TEMP");                                                                                 //create folder to hold temp MSA result
midway:	cout<<"Specify [ancthoring point threshold][gap amplification factor] seperate by space\n"                     //set position with a large proportion of identical AA as anchoring point
			<<"ep:90 3 represent [anchoring point threshold 90%][gap amplification factor 3]\n"
			<<"Input:";
		cin>>mancthld>>mamp;
		mtrial+=1;                                                                                               //control output path name
	   _Alignseqs(_seqsandinf._multiseqid,nbmatx,distmatx,_scoredb,seqsmatxuns,identical,_wlen,_gappnt,_cutoffper,_DorB,'Y',mancthld,mamp);
		score=_MSAsc(seqsmatxuns,_gappnt);
		MatxItoC(seqsmatxuns,dspmatx);                                                           
		string mwpath="c:\\PFSP\\TEMP\\Midway_Trial_"+inttostr(mtrial)+"_MSA_Preview_"+inttostr(mancthld)+"_"+inttostr(mamp)+".txt";  //midway path
		ofstream mwout(mwpath.c_str());
		mwout<<"This MSA result is for your preview only.Please delete this temp file in C:\\PFSP\\TEMP\n"
			 <<"Tips:Increase gap amplification factor to reduce gap.Use moderately lower anchor threshold to allow for more anchoring point\n" 
			 <<"Midway anchor threshold:"<<mancthld<<"  Midway gap amplification factor:"<<mamp<<"  Length:"<<dspmatx.getcnum()<<"  Score(Gonnet):"<<score<<endl;
		LandMarker(dspmatx.getcnum(),mwout);mwout<<endl;
		MileStone(dspmatx.getcnum(),mwout);mwout<<endl;
		mwout<<identical<<endl;
		mwout<<dspmatx;
		mwout<<identical<<endl;
		MileStone(dspmatx.getcnum(),mwout);mwout<<endl;
		LandMarker(dspmatx.getcnum(),mwout);mwout<<endl;
		mwout.close();
		TMPPATH=mwpath;                                                       //assign path file to global variable to justify multi-thread function
		DSPFILE(TMPPATH);
		cout<<"Are you satisfied with this MSA result?[Y/N]:";
		cin>>answer;
		if(answer=='N'||answer=='n')
			goto midway;
	}
	///////////////////////////////////////////////////////////////////////////Shift final MSA result to filter out obvious error
	int ancthld=INT_MIN,amp=INT_MIN;
	cout<<"Iterate final MSA result?[Y/N]:";
	cin>>answer;
	if(answer=='Y'||answer=='y')
       score=_MdfMSA(seqsmatxuns,dspmatx,identical,ancthld,amp,_gappnt);              //modify final MSA result,which might be less useful when midway iteration was enabled
	else if(score==(double)INT_MIN)
	   score=_MSAsc(seqsmatxuns,_gappnt);                                            //calculate MSA score,Gonnet scoring matrix was used
	
	seqsidmatx.fast_resize(_seqsandinf._seqnum,seqsmatxuns.getcnum());                                       //AAid seqs matrix,hold sorted AAid seqs after MSA
	for(int i=0;i<seqsmatxuns.getrnum();++i)                                                                 //resort seqs to the same as imported
	{
		for(int j=0;j<seqsmatxuns.getcnum();++j)
			seqsidmatx(seqsidorder[i],j)=seqsmatxuns(seqsmatxuns.getrnum()-i-1,j);
	}
    MatxItoC(seqsidmatx,dspmatx);
	cout<<"Output result..."<<endl;
	///////////////////////////////////////////////////////////////////////////////////////////////write log
	char time[10],date[10];
	_strtime(time);_strdate(date);
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<endl;                                           //process time was not computed since several steps need user option
	outlog<<"MSA "<<_seqsandinf._seqnum<<" sequences Distype:"<<((_EXorAP=='E'||_EXorAP=='e')?"Bayesian exact distance":"Aproximate distance")<<" Subtype:"<<_distype<<endl<<endl;
    ///////////////////////////////////////////////////////////////////////////////////////////////output detail result
	int width=(int)(log((double)dspmatx.getrnum())/log(10.0))+2;                                                 //format tool,dynamic ID marker length. 
	sarray<char> BLANK(width,' ');                                                                               //format tool
	int I;
	if(identical.size()==0)
	{
		identical.resize(dspmatx.getcnum());
		for(int j=0;j<dspmatx.getcnum();++j)
		{
			for(I=1;I<dspmatx.getrnum();++I)
			{
				if(dspmatx(I,j)!=dspmatx(0,j))
					break;
			}
			identical[j]=(I==dspmatx.getrnum()?'*':' ');
		}
	}
	MKDIR("C:\\PFSP\\MSA");
	string outpath,report,record;                                                                                //outpath hold common path
	outpath="c:\\PFSP\\MSA\\MSA_"+inttostr(dspmatx.getrnum())+"_"+string(&_EXorAP,sizeof(char))+"_"+string(&_distype,sizeof(char))+"_"+inttostr(_gappnt);
	report=outpath+".txt";                                                                                       //report hold path to detailed result
	ofstream outMSA(report.c_str());
	outMSA<<date<<" "<<time<<endl
		<<"Protein Functional Sites Prediction--MSA:(C++) Copyright 2009,Ruan Xiaoyang"<<endl
		<<"Multialign "<<_seqsandinf._seqnum<<" sequences  Distype:"<<((_EXorAP=='E'||_EXorAP=='e')?"Bayesian exact distance":"Aproximate distance")<<"  Subtype:"<<_distype
		<<"  Word length:"<<_wlen<<"  Gap penalty:"<<_gappnt<<"  Score(Gonnet):"<<score<<endl;
	if(mancthld!=INT_MIN && mamp!=INT_MIN)
		outMSA<<"Midway anchor threshold:"<<mancthld<<"  Midway gap amplification factor:"<<mamp<<"  ";
	if(ancthld!=INT_MIN && amp!=INT_MIN)
		outMSA<<"Anchor threshold:"<<ancthld<<"  Gap amplification factor:"<<amp;
	outMSA<<endl;
	if(_EXorAP=='E'||_EXorAP=='e')
	    outMSA<<"Bayesian step:"<<_step<<"  ";
	outMSA<<"Scoring database:DayhoffPAM1-400"<<endl<<endl;

	if(_seqsandinf._format=='F')
	{
		for(int i=0;i<_seqsandinf._seqnum;++i)
		{
			outMSA<<"ID ";
			outMSA.width(width);
			outMSA<<i<<" Gi_number:"<<_seqsandinf._infvec[i]._ginum<<" Accession:";
			outMSA<<_seqsandinf._infvec[i]._accession<<" Locus:"<<_seqsandinf._infvec[i]._locus<<endl;
		}
	}
	else if(_seqsandinf._format=='C')
	{
		for(int i=0;i<_seqsandinf._seqnum;++i)
		{
			outMSA<<"ID ";
			outMSA.width(width);
			outMSA<<i<<" "<<_seqsandinf._infvec[i]._gi<<endl;
		}
	}

    outMSA<<endl<<"Multiple sequence alignment result.'*'stands for identical"<<endl<<endl; 
	outMSA<<BLANK<<" ";LandMarker(dspmatx.getcnum(),outMSA);outMSA<<endl;
	outMSA<<BLANK<<" ";MileStone(dspmatx.getcnum(),outMSA);outMSA<<endl;
	outMSA<<BLANK<<" "<<identical<<endl;
	for(int i=0;i<dspmatx.getrnum();++i)
	{
		if(i!=0 && i!=(dspmatx.getrnum()-1) && i%50==0)
		{
			outMSA<<BLANK<<" ";LandMarker(dspmatx.getcnum(),outMSA);outMSA<<endl;
            outMSA<<BLANK<<" ";MileStone(dspmatx.getcnum(),outMSA);outMSA<<endl;
            outMSA<<BLANK<<" "<<identical<<endl;
		}
		outMSA.setf(ios::right);outMSA.width(width);
		outMSA<<i<<" "<<dspmatx[i]<<" "<<i<<endl;
	}
    outMSA<<BLANK<<" "<<identical<<endl;
	outMSA<<BLANK<<" ";MileStone(dspmatx.getcnum(),outMSA);outMSA<<endl;
	outMSA<<BLANK<<" ";LandMarker(dspmatx.getcnum(),outMSA);

	outMSA<<endl<<"Pairwise PAM distance matrix (PAM/100)"<<endl;
	for(int i=0;i<dspmatx.getrnum();++i)
	{
		outMSA<<"   ID ";
		outMSA.flags(ios::fixed);outMSA.width(4);outMSA.setf(ios::left);
		outMSA<<i;
	}
	outMSA<<endl<<distmatx<<endl;

	outMSA<<"Pairwise PAM distance variance matrix (PAM/100)"<<endl;
	for(int i=0;i<dspmatx.getrnum();++i)
	{
		outMSA<<"   ID ";
		outMSA.flags(ios::fixed);outMSA.width(4);outMSA.setf(ios::left);
		outMSA<<i;
	}
    outMSA<<endl<<varmatx<<endl;

	outMSA<<"Evolutionary tree.Values inside bracket are branch length (PAM/100)"<<endl;
	for(int i=0;i<nbmatx._nbsize;++i)
		outMSA<<nbmatx._marr1[i]<<"("<<nbmatx._len(i,0)<<")----("<<nbmatx._len(i,1)<<")"<<nbmatx._marr2[i]<<endl;

    record=outpath+"_record";                                                                                           //path to MSA record
	ofstream os(record.c_str());
	os<<dspmatx.getrnum()<<" "<<dspmatx.getcnum()<<endl;
	os<<identical<<endl;
	os<<dspmatx;
	outMSA.close();
	cout<<"MSA Result has been saved to "<<report<<endl;
	writeTStofile();
	TMPPATH=report;
	DSPFILE(TMPPATH);
	return;
}



