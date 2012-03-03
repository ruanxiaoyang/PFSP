vector<intmatrix> getwaypoint(darray<int> & _alignmatrix)
{	
	int i=0,j=0,k=0,fixsize;
	bool found;
	vector<intmatrix> allway;
	darray<int> coordinate(0,2,0);
	darray<int> temp(1,2,0);
	darray<int> addcoo(0);
	coordinate=_alignmatrix.maxwaypoint(0,0);
	while(i<coordinate.getrnum())
	{
		found=false;
	   if(coordinate(i,0)==0 || coordinate(i,1)==0)
       {		
			temp(0,0)=coordinate(i,0);
			temp(0,1)=coordinate(i,1);
			allway.push_back(temp);
	   }
		for(j=0;j<allway.size();++j)
		{
			if(allway[j](allway[j].getrnum()-1,0)==coordinate(i,0) && allway[j](allway[j].getrnum()-1,1)==coordinate(i,1))
			{
				found=true;
				break;
			}
		}
		if(found==true)
		{
			addcoo=_alignmatrix.maxwaypoint(coordinate(i,0)+1,coordinate(i,1)+1);
		    if(addcoo.getrnum()>0)
		    {
		         fixsize=allway.size();
		         for(j=0;j<fixsize;++j)
		         {
			        if(allway[j](allway[j].getrnum()-1,0)==coordinate(i,0) && allway[j](allway[j].getrnum()-1,1)==coordinate(i,1))
			        {
				    	allway[j].push_row(2);
                        allway[j](allway[j].getrnum()-1,0)=addcoo(0,0);
				        allway[j](allway[j].getrnum()-1,1)=addcoo(0,1);
				    	for(k=1;k<addcoo.getrnum();++k)
					    {
						   allway.push_back(allway[j]);
						   allway[allway.size()-1]((allway[allway.size()-1].getrnum())-1,0)=addcoo(k,0);
						   allway[allway.size()-1]((allway[allway.size()-1].getrnum())-1,1)=addcoo(k,1);
					    }
			         }
		         }
			coordinate.verticaladd(addcoo);
		    }
		}
		 i+=1;
	}
	return allway;
}

void createwaysdb(vector<intmatrix> & _aligndb,vector<vecintmatrix> & _waysdb)
{
	for(int i=0;i<_aligndb.size();++i)
		_waysdb.push_back(getwaypoint(_aligndb[i]));
}