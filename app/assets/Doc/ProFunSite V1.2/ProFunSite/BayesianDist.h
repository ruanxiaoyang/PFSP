//////////////////////////////////////////////////////////////
//  Bayesian Estimate of Exact Genetic Distance/Variance    //
//            Pankaj Agarwal and David J. States            //
//               J COMP BIO 3(1):1-17 1996                  //
//           Copyright Belongs to Ruan Xiaoyang             //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
//Parameter list instruction
//take 8 values
//"_aligntype" alignment type.'N/n'for Needleman global rule.'B/b'for BLAST local rule
//"_scoredb" standard Log10 based PAM score database(multiplied by 10).starts from PAM1,end at PAM400 or higher
//"_arrid1"/"_arrid2" protein id sequence data
//"_wlen" when "_aligntype" is 'B/b'(BLAST),this value is the word length of library (see RBLAST for detail)
//"_gappnt" gap penalty,this option is useful only when "_aligntype" is 'N/n'
//"_step" control how PAM scoring matrices were skipped.higher value means higher speed but higher variance also.the smallest value is 1,which means PAM scoring matrices will be applied one by one.  
//"_acculvl" accuracy level.a relative probability below 10^(-1*_acculvl) will be discarded,value higher than 4 gives almost same result.
//output 2 values
//"dist" hold distance value(PAM/100)
//"var" hold variance value(PAM/100)

void BayesianDist(const char & _aligntype,vector<darray<double>> &_scoredb,sarray<int> & _arrid1,sarray<int> & _arrid2,double & dist,double & var,const int & _wlen=3,const int & _gappnt=-3,const int & _step=3,const int & _acculvl=4)
{
	double score,maxscore;
	bool drop=false;
	int mark=0,acculvl=_acculvl*10;                                                            //since PAM score matrix was 10 time magnified
	double RO,ROtotal=0.0,sqE=0.0,prob;
	sarray<int> distarrtmp;
	sarray<double> ROarr,distarr;
	sarray<double> scorearr;
	darray<int> BKMATX,alnmatx;
	dist=0.0;var=0.0;
	if(_aligntype=='N' || _aligntype=='n')                                                     //if Needleman global alignment rule was specified,create global alignment matrix
	{
		NWalign(_arrid1,_arrid2,BKMATX);
        NWwp(BKMATX,_arrid1,_arrid2,Gonnet,_gappnt,alnmatx,score);                             //"Gonnet" score matrix.it was supposed to be the best score matrix for initial alignment.(science 1992 256:1443-1445)
	}
	else if(_aligntype=='B' || _aligntype=='b')                                                //if RBLAST Glocal alignment rule was specified,create glocal alignment matrix
		BLAST(_wlen,Gonnet,_arrid1,_arrid2,-1,alnmatx,95,'G');                                 //empirically fix word length to 3
	for(int i=0;i<_scoredb.size();i+=_step)
	{
		if(_aligntype=='N' || _aligntype=='n')  
			score=_Getscore(_scoredb[i],alnmatx,_arrid1.size(),_arrid2.size(),_arrid1,_arrid2,_gappnt);//see [RBLAST]-[_Getscore] 
		else if(_aligntype=='B' || _aligntype=='b')  
			score=_Getscore(_scoredb[i],alnmatx,_arrid1.size(),_arrid2.size(),_arrid1,_arrid2);
		if(drop==false && scorearr.size()>0)
		{
			if((score-scorearr[scorearr.size()-1])<0)
				mark+=1;
			else
				mark=0;
            if(mark==3)                                                                         //if 4 consecutive negative value was encountered,the value should be dropping
			{
				drop=true;
				maxscore=scorearr.smax();
			}
		}
		if(drop==true && ((maxscore>acculvl &&(maxscore-score)>acculvl)|| (maxscore<=acculvl && score<=0)))   
			break;	
		scorearr.pushback(score);
		distarrtmp.pushback(i);
	}
	if(drop==false)                                                                              //if no drop ever occurred,find the max score
		maxscore=scorearr.smax();
	if(maxscore-scorearr[scorearr.size()-1]<=acculvl)
		mark=scorearr.size()-1;
	else
	{
		for(int i=scorearr.size()-1;i>0;--i)
		{
			if(maxscore-scorearr[i]<=acculvl)
			{
				mark=i;
				break;
			}
		}
	}
    for(int i=0;i<scorearr.size();++i)
	{
		if(maxscore-scorearr[i]<=acculvl)                                                        //consider only those scores within accuracy threshold
		{
			RO=(double)pow(10.0,(double)(scorearr[i]-scorearr[mark])/10.0);                      
			ROarr.pushback(RO);
			distarr.pushback((double)(distarrtmp[i]+1));
			ROtotal+=RO;		                                                                 //total probability
		}
	}
    for(int i=0;i<ROarr.size();++i)
	{
		prob=ROarr[i]/ROtotal;                                                                   //relative probability
    	dist+=0.01*distarr[i]*prob;                                                              //Expectation E(X)
     	sqE+=pow(0.01*distarr[i],2.0)*prob;                                                      //Expectation of square value E(X^2)
	}
	var=sqE-pow(dist,2.0);                                                                       //variance is computed as E(X^2)-E(X)^2
}