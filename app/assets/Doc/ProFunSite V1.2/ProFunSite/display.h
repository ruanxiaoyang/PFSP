
//Display two aligned seqs."|"for identical."*"for similar(score>=0)
//Take 5 values."_arr1"/"_arr2" AA sequences."_scorematx"scoring matrix."_alnmatx"standard 0 or 1 matrix."_os"output stream.
template <class type>
void Dsptwo(sarray<char> & _arr1,sarray<char> & _arr2,darray<type> & _scorematx,darray<int> & _alnmatx,ostream & _os=cout)
{
	int rownum=_arr1.size(),colnum=_arr2.size(),previ=-1,prevj=-1,size,posi=-1,remain,gapi,gapj;
    size=Length(_alnmatx,rownum,colnum);                                          //compute length of sequence after alignment
	darray<char> multiseq(3,size,'a');
	for(int i=0;i<rownum;++i)
	{
		for(int j=prevj+1;j<colnum;++j)
		{
			if(_alnmatx(i,j)==1)
			{
				gapi=i-previ-1;
				gapj=j-prevj-1;
				if(gapi>0)
				{
					for(int g=1;g<=gapi;++g)
					{
						multiseq(0,posi+g)='-';
						multiseq(1,posi+g)=' ';
                        multiseq(2,posi+g)=_arr1[previ+g];
					}
					posi+=gapi;
				}
				if(gapj>0)
				{
					for(int g=1;g<=gapj;++g)
					{
						multiseq(0,posi+g)=_arr2[prevj+g];
						multiseq(1,posi+g)=' ';
						multiseq(2,posi+g)='-';
					}
					posi+=gapj;
				}
				posi+=1;
				previ=i;prevj=j;
				multiseq(0,posi)=_arr2[j];
				multiseq(1,posi)=(_arr2[j]==_arr1[i] ? '|' : (_scorematx(aaid(_arr2[j]),aaid(_arr1[i]))>0 ? '*' :' '));
				multiseq(2,posi)=_arr1[i];
				if(prevj==colnum-1)
				{
					if((remain=rownum-1-i)>0)
					{
						for(int g=1;g<=remain;++g)
						{
							multiseq(0,posi+g)='-';
							multiseq(1,posi+g)=' ';
							multiseq(2,posi+g)=_arr1[previ+g];
						}
					}
					LandMarker(size,_os);_os<<endl;
					MileStone(size,_os);_os<<endl;
					_os<<multiseq;
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
			multiseq(0,posi+g)=_arr2[prevj+g];
			multiseq(1,posi+g)=' ';
			multiseq(2,posi+g)='-';
		} 
	}
	LandMarker(size,_os);_os<<endl;
	MileStone(size,_os);_os<<endl;
	_os<<multiseq;
	return;
}
//Overloaded function of Dsptwo(sarray<char> & _arr1,sarray<char> & _arr2,...) to take AAid sequence (sarray<int> & _arrid1,sarray<int> & _arrid2....) as argument
template<class type>
void Dsptwo(sarray<int> & _arrid1,sarray<int> & _arrid2,darray<type> & _scorematx,darray<int> & _alnmatx,ostream & _os=cout)
{
	sarray<char> seq1(_arrid1.size()),seq2(_arrid2.size());
	for(int i=0;i<_arrid1.size();++i)
        seq1[i]=idaa(_arrid1[i]);
	for(int i=0;i<_arrid2.size();++i)
		seq2[i]=idaa(_arrid2[i]);
	Dsptwo(seq1,seq2,_scorematx,_alnmatx,_os);
}
