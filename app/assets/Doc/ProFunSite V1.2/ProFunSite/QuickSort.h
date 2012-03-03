//quick sort from min to max,<type> could be sarray or standard C array or vector
template <class type,class type2>
int partition(type & arr,int left,int right,type2 & temp)             //"temp" is indispensable for the compiler to know what type2 is
{
	int lo,hi;
    type2 pivot=arr[left],t;
    lo=left-1;
    hi=right+1;
    while(lo+1!=hi) 
	{
		if(arr[lo+1]<=pivot)
			lo++;
        else if(arr[hi-1]>pivot)
            hi--;
        else 
	    {
            t=arr[lo+1];
            arr[++lo]=arr[hi-1];
            arr[--hi]=t;
        }
    }
    arr[left]=arr[lo];
    arr[lo]=pivot;
    return lo;
}
template <class type>
void QuickSort(type & arr, int left,int right)
{
    int dp;
    if (left<right)
	{
		dp=partition(arr,left,right,arr[0]);
		QuickSort(arr,left,dp-1);
		QuickSort(arr,dp+1,right);
    }
}

