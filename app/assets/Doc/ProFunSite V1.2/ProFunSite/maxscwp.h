
//menu:[maxscwp]
//calculate standard 0 or 1 alignment matrix.Needleman global rule.
//take 5 values."_NWmatrix" Needleman alignment matrix."_scorematx" score matrix."_arr1"/"_arr2" AA sequence."_gappnt" gap penalty
//compute 2 value."alnmatx" standard 0 or 1 alignment matrix."alignscore" score of the final alignment.
template <class type>
void maxscwp(darray<int> & _NWmatrix,darray<type> & _scorematx,sarray<int> & _arrid1,sarray<int> & _arrid2,darray<int> & alnmatx,type & alignscore,const int & _gappnt=-3)
{
	int i=0,j=0,gap=0,mark,marki,markj;
	alignscore=(type)0;
	alnmatx.fast_resize(_arrid1.size(),_arrid2.size());
	darray<int> tempcoo;
	sarray<type> scorearr;
	sarray<int> maxscoreid,gaparr;
	while(i<_NWmatrix.getrnum() && j<_NWmatrix.getcnum())
	{
		marki=i;markj=j;
		tempcoo=_NWmatrix.maxwaypoint(i,j);                                                 //find the largest count in the way and store coordinate in tempcoo
		if(tempcoo.getrnum()==1)                                                            //if only one largest count found
		{
			i=tempcoo(0,0);j=tempcoo(0,1);
			alnmatx(i,j)=1;                                                                 //update "alnmatx"
			if(i==0)
			    alignscore+=_scorematx(_arrid1[i],_arrid2[j]);                            //if this is the start point,no gap penalty  
			else
				alignscore+=_scorematx(_arrid1[i],_arrid2[j])+(i-marki+j-markj)*_gappnt;  
			i+=1;j+=1;
			continue;
		}
		else if(tempcoo.getrnum()>1)                                                       //if more than one point have same largest count,choose the one with largest score
		{
			for(int k=0;k<tempcoo.getrnum();++k)        
			{
				gap=0;
				gap+=tempcoo(k,0)-i;gap+=tempcoo(k,1)-j;
				gaparr.pushback(gap);
				if(i!=0)
				   scorearr.pushback(_scorematx(_arrid1[tempcoo(k,0)],_arrid2[tempcoo(k,1)])+(type)gap*_gappnt);
				else
                   scorearr.pushback(_scorematx(_arrid1[tempcoo(k,0)],_arrid2[tempcoo(k,1)]));
			}
            alignscore+=scorearr.smax(maxscoreid);            
			scorearr.clear();
			if(maxscoreid.size()==1)
			{
				i=tempcoo(maxscoreid[0],0);j=tempcoo(maxscoreid[0],1);
				alnmatx(i,j)=1;
                i+=1;j+=1;
				gaparr.clear();
			    maxscoreid.clear();
				continue;
			}
			else                                                                           //if several points have same largest score,choose the one with smallest gap
			{
				mark=maxscoreid[0];
				for(int k=1;k<maxscoreid.size();++k)
				{
					if(i!=0)
					{
						if(gaparr[maxscoreid[k]]<gaparr[maxscoreid[k-1]])
						   mark=maxscoreid[k];
					}
					else if(i==0)
					{
						if(gaparr[maxscoreid[k]]>gaparr[maxscoreid[k-1]])
							mark=maxscoreid[k];
					}
				}
				i=tempcoo(mark,0);j=tempcoo(mark,1);
				alnmatx(i,j)=1;
				i+=1;j+=1;
				gaparr.clear();
			    maxscoreid.clear();
				continue;
			}
		}
	}
}





	


