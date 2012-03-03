//////////////////////////////////////////////////////////////
//                   Global Alignment                       //
//        SAUL B.NEEDLEMAN AND CHRISTIAN D.WUNSCH           //
//             J.MOL.BIOL.(1970) 48,443-453                 //
//           Copyright Belongs to Ruan Xiaoyang             //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
//This module performs global alignment between two amino acid sequences (nucleotide compatible)

template<class type>
inline void NWalign(sarray<type> & _arr1,sarray<type> & _arr2,darray<int> & BKMATX);

//Create alignment database
//Take 1 values."_multiseq" the multiple sequence database
//Update 2 value."_alignidvec" store pairing information of multiple sequences."_aligndb" Needleman alignment database for all sequence pairs
template<class type>
void NWALN(proformat & _seqsandinf,vector<darray<double>> &_scoredb,darray<type> & _scorematx,const int & _scmatxid,const char & _fixscmatx,const char & _subtype,const char & _DorB,const int & _gappnt,const int & _gma=1)
{
	MKDIR("C:\\PFSP\\NWALN");
	string cmnpath="C:\\PFSP\\NWALN\\NWALN_"+inttostr(_seqsandinf._seqnum)+"_"+string(&_fixscmatx,sizeof(char))+"_"+string(&_DorB,sizeof(char));
	string path;
	if(_fixscmatx=='Y'||_fixscmatx=='y')
	   path=cmnpath+"_"+inttostr(_scmatxid)+"_"+inttostr(_gappnt)+".txt";
	else
	   path=cmnpath+"_"+string(&_subtype,sizeof(char))+"_"+inttostr(_gappnt)+".txt";                    //name the path
	
	ofstream outNWA(path.c_str());
	char time[10],date[10];
	_strtime(time);_strdate(date);
    outNWA<<date<<" "<<time<<endl;
	outNWA<<"Protein Functional Sites Prediction--Needleman Global Alignment:(C++) Copyright 2009,Ruan Xiaoyang"<<endl;
	outNWA<<"Pairwise align "<<_seqsandinf._seqnum<<" sequences   Gap penalty:"<<_gappnt<<"  ";
	if(_fixscmatx=='Y' || _fixscmatx=='y')
		outNWA<<"Scoring matrix:"<<((_DorB=='D'|| _DorB=='d')?"PAM":"BLOSUM")<<_scmatxid<<endl<<endl;
	else
	{
		outNWA<<"Score database:DayhoffPAM1-400  Approximate distance:"<<_subtype;
		if(_subtype=='G' || _subtype=='g')
			outNWA<<"  Gamma:"<<_gma;
		outNWA<<endl<<endl;
	}
	clock_t start,end;                                                  //record process time for log file
	start=clock();
	int width=(int)(log((double)_seqsandinf._seqnum)/log(10.0))+2;      //set appropriate spacing
	double dist,var,smlt;                                               //distance,variance,similarity
	type alnscore;
	darray<int> BKMATX,alnmatx;                                         //background matrix,standard 0 or 1 alignment matrix
	cout<<"Needleman pairwise alignment in progress...";
	for(int i=0;i<_seqsandinf._seqnum;++i)                              //pairwise alignment
	{
		PctMarker(i,_seqsandinf._seqnum-1,2);
		for(int j=i+1;j<_seqsandinf._seqnum;++j)
		{
			if(_seqsandinf._format=='C')                               //when using custom format,only output _gi(stored all information inputed by user
			{
				    outNWA<<"ID";
					outNWA.width(width);
					outNWA<<j<<" "<<_seqsandinf._infvec[j]._gi<<endl;
					outNWA<<"ID";
					outNWA.width(width);
					outNWA<<i<<" "<<_seqsandinf._infvec[i]._gi<<endl;
			}
			else if(_seqsandinf._format=='F')
			{
					outNWA<<"ID";
					outNWA.width(width);
					outNWA<<j<<" Gi_number:"<<_seqsandinf._infvec[j]._ginum<<" Accession:"<<_seqsandinf._infvec[j]._accession<<" Locus:"<<_seqsandinf._infvec[j]._locus<<endl;
					outNWA<<"ID";
					outNWA.width(width);
					outNWA<<i<<" Gi_number:"<<_seqsandinf._infvec[i]._ginum<<" Accession:"<<_seqsandinf._infvec[i]._accession<<" Locus:"<<_seqsandinf._infvec[i]._locus<<endl;
			}
			if(_fixscmatx=='Y' || _fixscmatx=='y')
				NWalign(_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],BKMATX);
			else
			{
	            AprxDist(_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],dist,var,smlt,_subtype,BKMATX,_gma);
				outNWA<<"Estimated distance:PAM"<<(int)(dist*100)<<"  Variance:"<<var*100<<"  ";
				outNWA<<"Similarity:"<<smlt*100<<"%   ";
			}
			if(_fixscmatx=='Y' || _fixscmatx=='y')
			{
				NWwp(BKMATX,_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],_scorematx,_gappnt,alnmatx,alnscore);
				outNWA<<"Score(10*log[base10]):"<<alnscore<<endl;
			    Dsptwo(_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],_scorematx,alnmatx,outNWA);
			}
			else
			{
				NWwp(BKMATX,_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],_scoredb[(int)((dist*100)-1)],_gappnt,alnmatx,alnscore);
				outNWA<<"Score(10*log[base10]):"<<alnscore<<endl;
				Dsptwo(_seqsandinf._multiseqid[i],_seqsandinf._multiseqid[j],_scoredb[(int)((dist*100)-1)],alnmatx,outNWA);
			}
			outNWA<<endl;
		}
	}
	cout<<endl;
	outNWA.close();
	cout<<"NW Result has been saved to "<<path<<endl;
	end=clock();
	////////////////////////////////////////////////////////////////////////////////////////////write log
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<" Process time "<<end-start<<" minisec"<<endl;
	outlog<<"Global alignment of "<<_seqsandinf._seqnum<<" sequences"<<endl;
	if(_fixscmatx=='Y' || _fixscmatx=='y')
		outlog<<"Scoring matrix:"<<((_DorB=='D'|| _DorB=='d')?"PAM":"BLOSUM")<<_scmatxid<<"  ";
	else
		outlog<<"Score database:DayhoffPAM1-400  Approximate distance:"<<_subtype<<"  ";
	outlog<<"Gap penalty:"<<_gappnt<<endl<<endl;

	TMPPATH=path;                                               //display result
    DSPFILE(TMPPATH);
}

//Menu:[NWalign]
//Create Needleman alignment traceback matrix(same as the matrix generated by dynamic programming).It executes with high speed and low space occupation.It is faster than common dynamic programming
//Take 2 values."_arr1"/"_arr2" AA or AAid sequence (nucleotide compatible,<type> can be char or int)
//Return darray<int> Needleman alignment traceback matrix
template<class type>
void NWalign(sarray<type> & _arr1,sarray<type> & _arr2,darray<int> & BKMATX)
{
	int row=_arr1.size(),col=_arr2.size(),tmpsc,max,s;
	BKMATX.fast_resize(row+1,col+1,0);                          //BKMATX is a background matrix store the max match information.
	for(int j=col-1;j>=0;--j)
	{
		for(int i=row-1;i>=0;--i)
	    {		
			tmpsc=(_arr1[i]==_arr2[j]?1:0);
			s=tmpsc+BKMATX(i+1,j+1);                         
		    max=(BKMATX(i+1,j)>BKMATX(i,j+1)?BKMATX(i+1,j):BKMATX(i,j+1));
			BKMATX(i,j)=(s>max?s:max);
		} 
	}
	return;
}


//Menu:[NWalign]-[NWwp]
//This function is exclusively used to trace the waypoint of alignment matrix generated by NWalign function(Dynamic pogramming algorithm).
template<class type>
void NWwp(darray<int> & _BKMATX,sarray<int> & _arr1,sarray<int> & _arr2,darray<type> &_scmatx,const int & _gappnt,darray<int> & alnmatx,type & alnscore)
{
	int i=0,j=0,I,J,markI,markJ,gap,start,row=_arr1.size(),col=_arr2.size();
	alnscore=(type)0;
	type maxsc,tmpsc;
	alnmatx.fast_resize(row,col,0);
	while(i<row && j<col)
	{
		start=_BKMATX(i,j);
		maxsc=(type)INT_MIN;
		markI=i;markJ=j;
		if(_BKMATX(i+1,j+1)!=start)
		{
			J=j;
			while(J<col && _BKMATX(i,J)==start)
			{
			   if((tmpsc=_scmatx(_arr1[i],_arr2[J])+(J-j)*_gappnt)>maxsc)
			   {
				   markJ=J;
				   gap=markJ-j;
				   maxsc=tmpsc;
			   }
			   J+=1;
		   }
		   I=i+1;
		   while(I<row && _BKMATX(I,j)==start)
		   {
			   if((tmpsc=_scmatx(_arr1[I],_arr2[j])+(I-i)*_gappnt)>maxsc || (tmpsc==maxsc && (I-i)<gap))
			   {
				   markJ=j;
				   markI=I;
				   gap=I-i;
				   maxsc=tmpsc;
			   }
			   I+=1;
		   }
		}
		else
			maxsc=_scmatx(_arr1[markI],_arr2[markJ]);
	   if(markI==0 || markJ==0)
          alnscore+=_scmatx(_arr1[markI],_arr2[markJ]);
	   else
		   alnscore+=maxsc;
	   alnmatx(markI,markJ)=1;
	   i=markI+1;j=markJ+1;
	}
}



//NOT IN USE
//Dynamic programming.However,it is somewhat slower than the previous function.This might be caused by the long argument list in M()function
//Maybe a refinement on M() can improve the efficiency.But this can be done only by declaring NWmatx/arr1/arr2 as global variable. 
//template<class type>
//darray<int> NWalign(sarray<type> & _arr1,sarray<type> & _arr2)
//{
//	clock_t s,e;
//	s=clock();
//	darray<int> NWmatx(_arr1.size(),_arr2.size(),-1);
//	M(0,0,NWmatx,_arr1,_arr2);
//	e=clock();
//	cout<<"dynamic "<<e-s<<endl;
//	return NWmatx;
//}
//template<class type>
//int M(const int & _i,const int &_j,darray<int> &_NWmatx,sarray<type> & _arr1,sarray<type> & _arr2)
//{
//	if(_NWmatx(_i,_j)>=0)
//		return _NWmatx(_i,_j);
//	if(_i<_arr1.size()-1 && _j<_arr2.size()-1)
//	{
//		_NWmatx(_i,_j)=maximum(M(_i+1,_j+1,_NWmatx,_arr1,_arr2)+(_arr1[_i]==_arr2[_j]?1:0),M(_i+1,_j,_NWmatx,_arr1,_arr2),M(_i,_j+1,_NWmatx,_arr1,_arr2));
//		return _NWmatx(_i,_j);
//	}
//	if(_i==_arr1.size()-1 ||  _j==_arr2.size()-1)
//	{
//		_NWmatx(_i,_j)=(_arr1[_i]==_arr2[_j]?1:0);
//		return _NWmatx(_i,_j);
//	}
//}
//
//int maximum(const int & _a,const  int & _b,const int &_c)
//{
//	int max=(_a>_b? _a:_b);
//	return (max>_c?max:_c);
//}
//NWalign old version,slower than the one now used
//template<class type>
//darray<int> NWalign(sarray<type> & _arr1,sarray<type> & _arr2)
//{
//	int row=_arr1.size(),col=_arr2.size(),tmpsc,max;
//	darray<int> NWmatx(row,col),BKMATX(row,col);                                            //BKMATX is a background matrix store the max match information.It significantly speeds up the alignment
//	for(int j=col-1;j>=0;--j)
//	{
//		for(int i=row-1;i>=0;--i)
//	    {		
//			tmpsc=(_arr1[i]==_arr2[j]?1:0);
//			NWmatx(i,j)=tmpsc+((i!=row-1 && j!=col-1)?BKMATX(i+1,j+1):0);                   //only need to check BKMATX(i+1,j+1)to get the max value 
//			if(j==col-1)                                                                    //update BKMATX
//			{
//				if(i==row-1)
//					BKMATX(i,j)=tmpsc;
//				else
//					BKMATX(i,j)=(tmpsc>BKMATX(i+1,j) ? tmpsc:BKMATX(i+1,j));
//			}
//			else
//			{
//				if(i==row-1)
//					BKMATX(i,j)=(NWmatx(i,j)>BKMATX(i,j+1)?NWmatx(i,j):BKMATX(i,j+1));
//				else
//				{
//					max=(BKMATX(i+1,j)>BKMATX(i,j+1)?BKMATX(i+1,j):BKMATX(i,j+1));
//					BKMATX(i,j)=(NWmatx(i,j)>max?NWmatx(i,j):max);
//				}
//			}
//
//		} 
//	}
//	return NWmatx;
//}
//NWwp old version,slower than the one now used
//template<class type>
//void NWwp(darray<int> & _NWmatx,sarray<int> & _arr1,sarray<int> & _arr2,darray<type> &_scmatx,const int & _gappnt,darray<int> & alnmatx,type & alnscore)
//{
//	int i=0,j=0,marki,markj,row=_arr1.size(),col=_arr2.size(),_rindex_,_cindex_;
//	type rowmax,colmax,tmp;
//	alnmatx.fast_resize(row,col,0);
//	alnscore=(type)0;
//	darray<int> tempcoo;
//	while(i<row && j<col)
//	{
//		if(i==0 && j==0)
//		{
//			tempcoo=_NWmatx.maxwaypoint(i,j);
//			i=tempcoo(0,0);
//			j=tempcoo(0,1);
//			alnmatx(i,j)=1;
//			alnscore+=_scmatx(_arr1[i],_arr2[j]);
//			i+=1;
//			j+=1;
//		}
//		else
//		{
//			rowmax=(type)INT_MIN;
//			_NWmatx.rowmax(i,j,_cindex_);
//            for(int J=j;J<=_cindex_;++J)
//			{
//				tmp=_scmatx(_arr1[i],_arr2[J])+(type)((J-j)*_gappnt);
//				if(tmp>rowmax)
//				{
//					rowmax=tmp;
//					markj=J;
//				}
//			}
//			colmax=(type)INT_MIN;
//			_NWmatx.columnmax(j,i,_rindex_);
//			for(int I=i;I<=_rindex_;++I)
//			{
//				tmp=_scmatx(_arr1[I],_arr2[j])+(type)((I-i)*_gappnt);
//				if(tmp>colmax)
//				{
//					colmax=tmp;
//					marki=I;
//				}
//			}
//            if(rowmax>colmax || (rowmax==colmax && markj-j<=marki-i))
//			{
//				alnmatx(i,markj)=1;
//				alnscore+=rowmax;
//				i+=1;
//				j=markj+1;
//			}
//			else
//			{
//				alnmatx(marki,j)=1;
//				alnscore+=colmax;
//				i=marki+1;
//				j+=1;
//			}
//		}
//	}
//}