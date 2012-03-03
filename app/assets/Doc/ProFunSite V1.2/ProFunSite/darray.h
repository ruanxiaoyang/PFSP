//////////////////////////////////////////////////////////////
//                   Dynamic Matrix                         //
//           Copyright Belongs to Ruan Xiaoyang            //
//                    ruansun@163.com                       //
//////////////////////////////////////////////////////////////
using namespace std;

template <class type>
class darray
{
protected:int _row;
		int _column;
		type **dparray;
		type _defau;
public:darray();
	   darray(const type & defau);
	   darray(const darray<type> & temp);
	   darray(const int & row,const int & column);
	   darray(const int & row,const int & column,const type & defau);
	   void operator=(const darray<type> & temp);
	   darray<type> operator-(const darray<type> & temp);
	   type & operator()(const int & row,const int & column){return dparray[row][column];}
	   type rowmax(const int & row,const int & start,int & index);
	   type columnmax(const int & column,const int & start,int & index);
	   type columnmin(const int & column,const int & start,int & index);
	   type defau(){return _defau;}
	   int getrnum(){return _row;}
	   int getcnum(){return _column;}
	   int getcnum(const int & row);  
	   bool find(const int & row,const type & temp);
	   bool findc(const int & col,const type & temp);
	   sarray<type> getrow(const int & row);
	   sarray<type> getcol(const int & col);
	   sarray<type> operator[](const int & row);
	   darray<type> submatx(const int & i,const int & j,const int & k,const int & l);
	   darray<int> maxwaypoint(const int & row,const int & column);
	   darray<type> multiply(const darray<type> & temp);
	   darray<type> power(const int & power);
	   void verticalmerge(darray<type> & temp);
	   void verticaladd(darray<type> & temp);
	   void push_row(sarray<type> & temp);
	   void push_row(const int & column);
	   void push_to_row(const int & row,sarray<type> & temp);
	   void push_to_row(const int & row,const type & temp);
	   void poprow(const int & num_row);
	   void resize(const int & row,const int & column);
	   void writecolnum(const int & row);
	   void fast_resize(const int & row,const int & column);
	   void fast_resize(const int & row, const int &column,const type & defau);
	   void traverse();
	   void freem(){for(int i=0;i<_row;++i){free(dparray[i]);dparray[i]=NULL;}free(dparray);dparray=NULL;_column=0;_row=0;}
	   void clear();
	   friend ostream & operator<< <type>(ostream & os,darray<type> & temp); 
	   void record(ostream & os);
	   ~darray();
};

template<class type>
darray<type>::~darray()
{
	this->freem();
}
template <class type>
void darray<type>::clear()
{
	for(int i=0;i<_row;++i)
	{
		free(dparray[i]);
		dparray[i]=NULL;
	}
	dparray=(type**)realloc(dparray,0);
	_column=0;
	_row=0;
}
template<class type>
void darray<type>::writecolnum(const int & row)
{
	_column=this->getcnum(row);
}
template<class type>
darray<type> darray<type>::operator-(const darray<type> & temp)
{
	if(_row!=temp._row || _column!=temp._column)
	{
		cout<<"Minus operation failed!Return *this\n";
		return *this;
	}
	darray<type> result(_row,_column);
	for(int i=0;i<_row;++i)
	{
		for(int j=0;j<_column;++j)
			result(i,j)=dparray[i][j]-temp.dparray[i][j];
	}
	return result;
}
template<class type>
void darray<type>::poprow(const int & num_row)
{
	if(num_row>_row)
		return;
	_row-=num_row;
	for(int i=_row;i<_row+num_row;++i)
	{
		free(dparray[i]);
		dparray[i]=NULL;
	}
	dparray=(type**)realloc(dparray,sizeof(type*)*_row);
}
template<class type>
darray<type>::darray()
{
	_row=0;
	_column=0;
    _defau=(type)0;
    dparray=(type**)malloc(0);
}
template <class type>
darray<type>::darray(const type & defau)
{
	_row=0;
	_column=0;
	_defau=defau;
    dparray=(type**)malloc(0);
}
template<class type>
darray<type>::darray(const darray<type> &temp)
{
	if(this==&temp)
		return;
	_row=temp._row;
	_column=temp._column;
	_defau=temp._defau;
	dparray=(type**)malloc(sizeof(type*)*_row);
	for(int i=0;i<_row;++i)
		dparray[i]=(type*)malloc(sizeof(type)*_column);
	for(int i=0;i<_row;++i)
		for(int j=0;j<_msize((type*)temp.dparray[i])/sizeof(type);++j)
			dparray[i][j]=(temp.dparray)[i][j];
}
template<class type>
darray<type>::darray(const int & row,const int & column,const type & defau)
{
	_row=row;
	_column=column;
	_defau=defau;
    dparray=(type**)malloc(sizeof(type*)*_row);
	for(int i=0;i<_row;++i)
	   dparray[i]=(type*)malloc(sizeof(type)*_column);
	for(int i=0;i<_row;++i)
	    for(int j=0;j<_column;++j)
			dparray[i][j]=_defau;
}
template<class type>
darray<type>::darray(const int & row,const int & column)
{
	_row=row;
	_column=column;
	_defau=(type)0;
    dparray=(type**)malloc(sizeof(type*)*_row);
	for(int i=0;i<_row;++i)
	   dparray[i]=(type*)malloc(sizeof(type)*_column);
}
template<class type>
bool darray<type>::find(const int & row,const type & temp)
{
	for(int j=0;j<_msize((type*)dparray[row])/sizeof(type);++j)
	{
		if(dparray[row][j]==temp)
			return true;
	}
	return false;
}
template<class type>
bool darray<type>::findc(const int & col,const type & temp)
{
	for(int i=0;i<_row;++i)
	{
		if(dparray[i][col]==temp)
			return true;
	}
	return false;
}
template<class type>
darray<type> darray<type>::power(const int &power)
{
	darray<type> result,temp(*this);
	for(int i=0;i<power;++i)
	{
		result=temp.multiply(*this);
		temp=result;
	}
	return result;
}

template<class type>
darray<type> darray<type>::multiply(const darray<type> & temp)
{
    if(_column!=temp._row)
	{
		cerr<<"Row & Column number mismatch!\n";
     	return *this;
	}
	type p=(type)0;
	darray<type> result(_row,temp._column);
	for(int i=0;i<_row;++i)
	{
		for(int j=0;j<temp._column;++j)
		{
			for(int k=0;k<_col;++k)
				p+=dparray[i][k]*temp.dparray[k][j];
			result(i,j)=p;
			p=(type)0;
		}
	}
    return result;
}

template<class type>
void darray<type>::push_to_row(const int & row,sarray<type> & temp)
{
	int J=this->getcnum(row),j=0;
	int column=J+temp.size();
	dparray[row]=(type*)realloc(dparray[row],sizeof(type)*column);
	while(J<column)
	{
		dparray[row][J]=temp[j];
		J+=1;j+=1;
	}
}
template<class type>
void darray<type>::push_to_row(const int & row,const type & temp)
{
	int column=this->getcnum(row)+1;
	dparray[row]=(type*)realloc(dparray[row],sizeof(type)*column);
	dparray[row][column-1]=temp;
}
template <class type>
int darray<type>::getcnum(const int &row)
{
     return _msize((type*)dparray[row])/sizeof(type);
}
template<class type>
void darray<type>::push_row(const int & column)
{
	_row+=1;
	if(column>_column)
	   _column=column;
	dparray=(type **)realloc(dparray,sizeof(type *)*_row);
    dparray[_row-1]=(type*)malloc(sizeof(type)*column);
}
template <class type>
void darray<type>::push_row(sarray<type> &temp)
{
	if(temp.size()==0)
		return;
	this->push_row(temp.size());
	for(int i=0;i<temp.size();++i)
		dparray[_row-1][i]=temp[i];
}
template <class type>
void darray<type>::verticaladd(darray<type> &temp)
{
	if(_column!=temp._column)
		return;
	for(int i=0;i<temp._row;++i)
	{
		this->push_row(_column);
        for(int l=0;l<_column;++l)
			dparray[_row-1][l]=temp.dparray[i][l];
	}
}

template <class type>
void darray<type>::verticalmerge(darray<type> &temp)
{
	if(_column!=temp._column)
		return;
	int mark=0,r=0,j=0;
	for(int i=0;i<temp._row;++i)
	{      
		j=0;
        while(j<_column && _row>0)
		{
			if(temp(i,j)==dparray[r][j])
			{ 
				mark+=1;j+=1;
				if(mark==_column)
					goto nextcom;
			    continue;
			}
			if(temp(i,j)!=dparray[r][j])
			{ 	
			       r+=1;
				   j=0;
				   mark=0;
				if(r==_row)
				goto addrow;
			}
		}
addrow:		this->push_row(_column);
		    for(int l=0;l<_column;++l)
			dparray[_row-1][l]=temp(i,l);
nextcom:	mark=0;	r=0;			
	}
}
template <class type>
darray<int> darray<type>::maxwaypoint(const int &row, const int &col)
{
	int colid,rowid;
	type rowmax,colmax;
	rowmax=this->rowmax(row,col,colid);
	if(row<_row-1)
		colmax=this->columnmax(col,row+1,rowid);
	else
		colmax=(type)INT_MIN;
	darray<int> wpmatx(1,2);
    if(rowmax>colmax)
	{
		wpmatx(wpmatx.getrnum()-1,0)=row;
		wpmatx(wpmatx.getrnum()-1,1)=colid;
	    for(int j=colid+1;j<_column;++j)
	    {
		    if(dparray[row][j]==rowmax)
		    {
				wpmatx.push_row(2);
				wpmatx(wpmatx.getrnum()-1,0)=row;
				wpmatx(wpmatx.getrnum()-1,1)=j;
		    }
		}
		return wpmatx;
	}
	if(colmax>rowmax)
	{
		wpmatx(wpmatx.getrnum()-1,0)=rowid;
		wpmatx(wpmatx.getrnum()-1,1)=col;
		for(int i=rowid+1;i<_row;++i)
	    {
		    if(dparray[i][col]==colmax)
		    {
				wpmatx.push_row(2);
				wpmatx(wpmatx.getrnum()-1,0)=i;
				wpmatx(wpmatx.getrnum()-1,1)=col;
		    }
	    }
		return wpmatx;
	}
	if(rowmax==colmax)
	{
        wpmatx(wpmatx.getrnum()-1,0)=row;
		wpmatx(wpmatx.getrnum()-1,1)=colid;
		if(row<_row-1)
		{
		   wpmatx.push_row(2);
           wpmatx(wpmatx.getrnum()-1,0)=rowid;
		   wpmatx(wpmatx.getrnum()-1,1)=col;
		}
		for(int j=colid+1;j<_column;++j)
		{
			if(dparray[row][j]==rowmax)
			{
				wpmatx.push_row(2);
                wpmatx(wpmatx.getrnum()-1,0)=row;
		        wpmatx(wpmatx.getrnum()-1,1)=j;
			}
		}
		if(row<_row-1)
		{
			for(int i=rowid+1;i<_row;++i)
		    {
				if(dparray[i][col]==colmax)
				{
					wpmatx.push_row(2);
					wpmatx(wpmatx.getrnum()-1,0)=i;
					wpmatx(wpmatx.getrnum()-1,1)=col;
				}
		    }
		}
        return wpmatx;
	}
}

template <class type>
type darray<type>::rowmax(const int &row,const int &start,int & index)
{
	type max=dparray[row][start];
	index=start;
	for(int j=start+1;j<_column;++j)
	{
		if(max<dparray[row][j])
		{
			max=dparray[row][j];
            index=j;
		}
	}
    return max;
}
template <class type>
type darray<type>::columnmax(const int &column,const int &start,int & index)
{
	type max=dparray[start][column];
	index=start;
	for(int i=start+1;i<_row;++i)
	{
		if(max<dparray[i][column])
		{
			max=dparray[i][column];
			index=i;
		}
	}
	return max;
}
template<class type>
type darray<type>::columnmin(const int & column,const int & start,int & index)
{
	type min=dparray[start][column];
	index=start;
	for(int i=start+1;i<_row;++i)
	{
		if(min>dparray[i][column])
		{
			min=dparray[i][column];
			index=i;
		}
	}
	return min;
}

template <class type>
darray<type> darray<type>::submatx(const int &i, const int &j, const int &k, const int &l)
{
	if(i<0 || i>_row-1 || j<0 || j>_column-1 || k<=i || k>_row-1 || l<=j || l>_column-1)
	{
		cout<<"Coordinate error,Submatx operation failed!\n";
		return *this;
	}
	int row=k-i+1;
	int column=l-j+1;
	darray<type> matx(row,column);
	for(int m=0;m<row;++m)
		for(int n=0;n<column;++n)
			matx(m,n)=dparray[i+m][j+n];
	return matx;
}
template <class type>
sarray<type> darray<type>::getrow(const int &row)
{
	int column=_msize((type*)dparray[row])/sizeof(type);
	sarray<type> temp(column);
	for(int j=0;j<column;++j)
		temp[j]=dparray[row][j];
	return temp;
}
template <class type>
sarray<type> darray<type>::operator[](const int &row)
{
	int column=_msize((type*)dparray[row])/sizeof(type);
	sarray<type> temp(column);
	for(int j=0;j<column;++j)
		temp[j]=dparray[row][j];
	return temp;
}
template <class type>
sarray<type> darray<type>::getcol(const int &col)
{
	sarray<type> temp(_row);
	for(int i=0;i<_row;++i)
		temp[i]=dparray[i][col];
	return temp;
}

template<class type>
void darray<type>::resize(const int & row,const int & column)
{
	darray<type> temp(*this);
	this->freem();
	_row=row;
	_column=column;
	dparray=(type**)malloc(sizeof(type*)*_row);
	for(int i=0;i<_row;++i)
		dparray[i]=(type*)malloc(sizeof(type)*_column);
	for(int i=0;i<((_row>temp._row) ? temp._row:_row);++i)
		for(int j=0;j<((_column>temp._column) ? temp._column:_column);++j)
			dparray[i][j]=temp(i,j);
    if(_row>temp._row)
	{
		for(int i=temp._row;i<_row;++i)
		for(int j=0;j<_column;++j)
			dparray[i][j]=_defau;
	}
	if(_column>temp._column) 
	{
		for(int i=0;i<((_row>temp._row) ? temp._row:_row);++i)
			for(int j=temp._column;j<_column;++j)
            dparray[i][j]=_defau;
	}
}
template<class type>
void darray<type>::fast_resize(const int &row, const int &column)
{
	this->freem();
	_row=row;
	_column=column;
	dparray=(type**)malloc(sizeof(type*)*_row);
	for(int i=0;i<_row;++i)
		dparray[i]=(type*)malloc(sizeof(type)*_column);
}
template<class type>
void darray<type>::fast_resize(const int & row, const int &column,const type & defau)
{
	this->freem();
	_row=row;
	_column=column;
	_defau=defau;
	dparray=(type**)malloc(sizeof(type*)*_row);
	for(int i=0;i<_row;++i)
	{
		dparray[i]=(type*)malloc(sizeof(type)*_column);
		for(int j=0;j<_column;++j)
			dparray[i][j]=_defau;
	}
}

template<class type>
void darray<type>::operator =(const darray<type> & temp)
{
	if(this==&temp)
		return;
	this->freem();
	_row=temp._row;
	_column=temp._column;
	_defau=temp._defau;
	int column;
    dparray=(type**)malloc(sizeof(type*)*_row);
	for(int i=0;i<_row;++i)
	{
		column=_msize((type*)(temp.dparray[i]))/sizeof(type);
		dparray[i]=(type*)malloc(sizeof(type)*column);
	}
	for(int i=0;i<_row;++i)
	    for(int j=0;j<(_msize((type*)(temp.dparray[i]))/sizeof(type));++j)
			dparray[i][j]=(temp.dparray)[i][j];
}
template<class type>
void darray<type>::record(ostream &os)
{
	if(_row==0)
	{
		os<<"Empty array!\n";
		return;
	}
    for(int i=0;i<_row;++i)
	{
		for(int j=0;j<this->getcnum(i);++j)
			os<<dparray[i][j]<<" ";
		os<<endl;
	}
}
template<class type>
ostream & operator<<(ostream & os,darray<type> & temp)
{
	if(temp._row==0)
	{	
		os<<"Empty array!\n";
	    return os;
	}
	int row=(_msize((type**)temp.dparray))/(sizeof(type*));
	if(sizeof(type)>4)
	{
	    os.flags(ios::fixed);os.precision(3);
	    for(int i=0;i<row;++i)
	    {
		   int column=(_msize((type*)temp.dparray[i]))/(sizeof(type));
		   for(int j=0;j<column;++j)
		   {
		    	os.width(9);
			    os<<(temp.dparray)[i][j]<<" ";
	       }
		   os<<endl;
	    }
	}
	else if(sizeof(type)==4)
	{
		for(int i=0;i<row;++i)
	    {
		   int column=(_msize((type*)temp.dparray[i]))/(sizeof(type));
		   for(int j=0;j<column;++j)
			   os<<(temp.dparray)[i][j]<<" ";
		   os<<endl;
		}
	}
	else
	{
		for(int i=0;i<row;++i)
	    {
		   int column=(_msize((type*)temp.dparray[i]))/(sizeof(type));
		   for(int j=0;j<column;++j)
			   os<<(temp.dparray)[i][j];
		   os<<endl;
		}
	}
	return os;
}
template<class type>
void darray<type>::traverse()
{
	darray<type> temp(_column,_row);
	for(int i=0;i<temp._row;++i)
		for(int j=0;j<temp._column;++j)
			temp.dparray[i][j]=dparray[j][i];
	*this=temp;
}