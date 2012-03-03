//////////////////////////////////////////////////////////////
//                 Neighbor-Joining Method                  //
//             Naruya Saitou and Masatoshi Nei              //
//             Mol.Biol.Evol.4(4):406-425.1987              //
//           Copyright Belongs to Ruan Xiaoyang             //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
//This module output unrooted genetic tree of multiple seqs,which guide the MSA process 
//For each run,a log file was saved at c:\ProFunSite.log indicating the tree information

//structure hold topology of multiple alignment seqs
struct nbstr
{
	darray<int> _marr1;                                      //left neighbour id
	darray<int> _marr2;                                      //right neighbour id
	darray<double> _len;                                     //store length of left and right neighbout to intermidiate node
	darray<int> _lkbk;                                       //look back previously aligned seqs.this can speed up the alignment.
	int _nbsize;
};
//Calculate genetic tree
//Take 1 value."_distmatx" distance matrix.
//Compute 1 value."nbmatx" neighbour matrix,which specify the topology of genetic tree and related distance
void NJtree(darray<double> & _distmatx,nbstr & nbmatx)
{    clock_t start,end;
     start=clock();
	 nbmatx._nbsize=0;                                                                                     //initialize to 0
	 darray<double> distmatx(_distmatx),distmatxnew,rough_len(_distmatx.getrnum()-2,2);                    //copy "_distmatx" to "distmatx"
     sarray<double> Dkarr(distmatx.getrnum()),Dkarrnew;                                                    //"Dkarr" array hold Dk value(Dk:sum over kth row of distance matrix devided by row-2)
	 double temp=0;
	 for(int i=0;i<distmatx.getrnum();++i)                                                                 //initiate Dkarr(in the following joining process very little computation is needed to update Dkarr)
	 {
		 temp=0;
		 for(int j=0;j<distmatx.getcnum();++j)
		 {
			 if(i!=j)
				 temp+=distmatx(i,j);
		 }
		 Dkarr[i]=temp/(distmatx.getrnum()-2);
	 }
	 darray<int> idarray(distmatx.getrnum(),1,0),idarraynew;                                               //store booking information of id
	 for(int i=0;i<distmatx.getrnum();++i)
         idarray(i,0)=i;
	 int rownum,marki,markj,rowiter,coliter;
	 double Q,Qmin=DBL_MAX,Lix,Ljx;
	 cout<<"Finding neighbour...";
	 while(distmatx.getrnum()>3)                                                                          //stop until 3 neighbours left
	 {
		 rownum=distmatx.getrnum();
		 PctMarker(_distmatx.getrnum()-rownum,_distmatx.getrnum()-4,2);                                   //progressbar
		 Qmin=DBL_MAX;
		 for(int i=0;i<Dkarr.size();++i)
		 {
			 for(int j=i+1;j<Dkarr.size();++j)
			 {
				 Q=distmatx(i,j)-Dkarr[i]-Dkarr[j];                                                       //the nature of finding a neighbour is to find the minimum Q,which maximize the common path linking i,j and other nodes
				 if(Qmin>Q)
				 {
					 Qmin=Q;
					 marki=i;
					 markj=j;
				 }
			 }
		 }
		Lix=(distmatx(marki,markj)+Dkarr[marki]-Dkarr[markj])/2.0;                                                   //distance of terminal i to node x
		Ljx=(distmatx(marki,markj)+Dkarr[markj]-Dkarr[marki])/2.0;                                                   //distance of terminal j to node x
////////////////////////////////////////////////////////////////construct neighbour 
		nbmatx._marr1.push_row(idarray[marki]);
		nbmatx._marr2.push_row(idarray[markj]);
		nbmatx._nbsize+=1;
////////////////////////////////////////////////////////////// 
//calculate rough length between sequences.since the "terminal" may be comprised of more sub terminals,so true length can not be directly calculated in this step 
		rough_len(_distmatx.getrnum()-rownum,0)=Lix;
		rough_len(_distmatx.getrnum()-rownum,1)=Ljx;
///////////////////////////////////////////////////////////////update distance matrix
        distmatxnew.fast_resize(rownum-1,rownum-1);
		distmatxnew(0,0)=0.0;
		coliter=1;
        for(int k=0;k<rownum;++k)
		{
			if(k!=marki && k!=markj)
			{
				distmatxnew(coliter,0)=distmatxnew(0,coliter)=(distmatx(marki,k)+distmatx(markj,k)-distmatx(marki,markj))/2.0;   //distance between newly formed node and other nodes
				distmatxnew(coliter,coliter)=0.0;
				coliter+=1;
			}
		}
		rowiter=1;
	    for(int i=0;i<rownum;++i)
		{
			if(i!=marki && i!=markj)
			{
				coliter=1;
				for(int j=0;j<rownum;++j)
				{
					if(j!=marki && j!=markj)
					{
						distmatxnew(rowiter,coliter)=distmatx(i,j);
						coliter+=1;
					}
				}
				rowiter+=1;
			}
		}
		temp=0.0;
		int iter=1;
        Dkarrnew.resize(rownum-1);                                                                            //update Dk arrary
		for(int j=0;j<distmatxnew.getrnum();++j)
            temp+=distmatxnew(0,j);
		Dkarrnew[0]=temp/(distmatxnew.getrnum()-2);
		for(int i=0;i<Dkarr.size();++i)
		{
			if(i!=marki && i!=markj)
			{
				Dkarrnew[iter]=(Dkarr[i]*(rownum-2)-distmatx(marki,i)-distmatx(markj,i)+distmatxnew(iter,0))/(double)(rownum-3);   //we need not to recompute Dk from new distance matrix,but rather can obtain it in an easy way
				iter+=1;
			}
		}
		distmatx=distmatxnew;
		Dkarr=Dkarrnew;
////////////////////////////////////////////////////////////////update sequence IDs
		idarraynew.push_row(idarray[marki]);
		idarraynew.push_to_row(idarraynew.getrnum()-1,idarray[markj]);
		for(int i=0;i<idarray.getrnum();++i)
		{
			if(i!=marki && i!=markj)
				idarraynew.push_row(idarray[i]);
		}
		idarray=idarraynew;
		idarraynew.clear();
	}                                                         //while loop stopped here
////////////////////////////////////////////////////////////////process ID information of the final three nodes
	nbmatx._marr1.push_row(idarray[1]);
	nbmatx._marr2.push_row(idarray[2]);
	nbmatx._nbsize+=1;
	nbmatx._marr1.push_row(idarray[0]);
	nbmatx._marr2.push_row(idarray[1]);
	nbmatx._marr2.push_to_row(nbmatx._nbsize,idarray[2]);
	nbmatx._nbsize+=1;
////////////////////////////////////////////////////////////////process rough length of the final two nodes
	rough_len(rough_len.getrnum()-1,0)=(distmatx(1,2)+distmatx(1,0)-distmatx(2,0))/2.0;
    rough_len(rough_len.getrnum()-1,1)=(distmatx(1,2)+distmatx(2,0)-distmatx(1,0))/2.0;
////////////////////////////////////////////////////////////////process real length from rough length and construct look back information
	nbmatx._len.fast_resize(nbmatx._nbsize,2);                //the size of _len and _lkbk is same as _marr1/_marr2 
	nbmatx._lkbk.fast_resize(nbmatx._nbsize,2);
	cout<<endl<<"Calculate length...";
	int r,r1,r2,colnum1,colnum2;
	for(r=0;r<nbmatx._nbsize-1;++r)
	{
		PctMarker(r,nbmatx._nbsize-2,2);
		colnum1=nbmatx._marr1.getcnum(r);
		colnum2=nbmatx._marr2.getcnum(r);
		if(colnum1==1 && colnum2==1)                          //if both _marr1 and _marr2 are non-aligned singal seq
		{
			nbmatx._len(r,0)=rough_len(r,0);                  //rough length is equal to real length(since there is no sub node) 
			nbmatx._len(r,1)=rough_len(r,1);
			nbmatx._lkbk(r,0)=-1;                             //no need to find previous record 
			nbmatx._lkbk(r,1)=-1;
		}
		else if(colnum1>1 && colnum2==1)                      //if _marr1 is aligned and _marr2 is non-aligned single seq 
		{
			for(r1=0;r1<r;++r1)
			{
				if(nbmatx._marr1[r]==(nbmatx._marr1[r1]+nbmatx._marr2[r1]))           //find previous record for match
					break;
			}
			nbmatx._len(r,0)=rough_len(r,0)-(rough_len(r1,0)+rough_len(r1,1))/2.0;    //process rough length into real length
			nbmatx._len(r,1)=rough_len(r,1);
			nbmatx._lkbk(r,0)=r1;                             //this means we need to use NO.r1 record (aligned seqs) to complete the current alignment
			nbmatx._lkbk(r,1)=-1;
		}
		else if(colnum1==1 && colnum2>1)
		{
			for(r2=0;r2<r;++r2)
			{
				if(nbmatx._marr2[r]==(nbmatx._marr1[r2]+nbmatx._marr2[r2]))
					break;
			}
			nbmatx._len(r,0)=rough_len(r,0);
			nbmatx._len(r,1)=rough_len(r,1)-(rough_len(r2,0)+rough_len(r2,1))/2.0;
			nbmatx._lkbk(r,0)=-1;
			nbmatx._lkbk(r,1)=r2;
		}
		else if(colnum1>1 && colnum2>1)
		{
			for(r1=0;r1<r;++r1)
			{
				if(nbmatx._marr1[r]==(nbmatx._marr1[r1]+nbmatx._marr2[r1]))
					break;
			}
			for(r2=0;r2<r;++r2)
			{
				if(nbmatx._marr2[r]==(nbmatx._marr1[r2]+nbmatx._marr2[r2]))
					break;
			}
			nbmatx._len(r,0)=rough_len(r,0)-(rough_len(r1,0)+rough_len(r1,1))/2.0;
			nbmatx._len(r,1)=rough_len(r,1)-(rough_len(r2,0)+rough_len(r2,1))/2.0;
			nbmatx._lkbk(r,0)=r1;
			nbmatx._lkbk(r,1)=r2;
		}
	}
	nbmatx._len(nbmatx._nbsize-1,0)=(distmatx(0,1)+distmatx(0,2)-(_distmatx.getrnum()>3? rough_len(rough_len.getrnum()-2,0)+rough_len(rough_len.getrnum()-2,1):0.0)-nbmatx._len(nbmatx._nbsize-2,0)-nbmatx._len(nbmatx._nbsize-2,1))/2.0;//process the final length
	nbmatx._len(nbmatx._nbsize-1,1)=0;
	nbmatx._lkbk(nbmatx._nbsize-1,0)=nbmatx._nbsize-3;             //the final neighbour was composed of the last two records 
	nbmatx._lkbk(nbmatx._nbsize-1,1)=nbmatx._nbsize-2;
	end=clock();
/////////////////////////////////////////////////////////////////////write log file
	char time[10],date[10];
	_strtime(time);_strdate(date);
	ofstream outlog("c:\\ProFunSite.log",ios_base::app);
    outlog<<date<<" "<<time<<" Process time "<<end-start<<" minisec"<<endl;
	outlog<<"NJtree processed "<<_distmatx.getrnum()<<" Sequences.Construced genetic tree with "<<nbmatx._nbsize<<" neighbours"<<endl<<endl;
}