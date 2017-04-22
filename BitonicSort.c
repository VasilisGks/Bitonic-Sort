#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <math.h>
#include <omp.h>
#define MAX 2

typedef struct Ar {
    int low;
    int high;
	int dir;
	int thd;
} element;
pthread_mutex_t mute,mute1;
struct timeval startwtime, endwtime;
double seq_time;
//*********************************************
volatile  cout1=0,cout2=0;
//pthread_mutex_t mute=PTHREAD_MUTEX_INITIALIZER;

int flag=1;
int q;           //thread max number
int N;          // data array size
int *a;         // data array to be sorted
//***********************************************
const int ASCENDING  = 1;
const int DESCENDING = 0;


void init(void);
void print(void);
void sort(void);
void test(void);
inline void exchange(int i, int j);
void compare(int i, int j, int dir);
void bitonicMerge(int lo, int cnt, int dir);
void recBitonicSort(int lo, int cnt, int dir);
void impBitonicSort(void);
//new
int CompAsc(const void *x,const void *y);
int CompDesc(const void *x,const void *y);
void impBitonicSort(void);
void *p_recBitonicSort(void *arr);
void *p_bitonicMerge(void *ar);
void p_sort(void);

/** the main program **/ 
int main(int argc, char **argv) {

  if (argc != 3) {
    printf("Usage: %s q\n  where n=2^q is problem size (power of two)\n", 
	   argv[0]);
    exit(1);
  }

  N = 1<<atoi(argv[1]);
  a = (int *) malloc(N * sizeof(int));
  q= 1<<atoi(argv[2]);                     //q is number of threads!!!
  
  init();

  gettimeofday (&startwtime, NULL);
  impBitonicSort();	
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Imperative wall clock time = %f\n", seq_time);

  test();

  init();
  gettimeofday (&startwtime, NULL);
  sort();
  gettimeofday (&endwtime, NULL);

  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Recursive wall clock time = %f\n", seq_time);

  test();
  //PARALLEL ************************************
  init();
  gettimeofday (&startwtime, NULL);	
	p_sort();
	 gettimeofday (&endwtime, NULL);
	 
  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Parallel wall clock time = %f\n", seq_time);
	test();
  // print();
  
  //qsort()
  
   init();
  gettimeofday (&startwtime, NULL);	
	qsort(a,N,sizeof(int),CompAsc);
	 gettimeofday (&endwtime, NULL);
	 
  seq_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6
		      + endwtime.tv_sec - startwtime.tv_sec);

  printf("Quicksort wall clock time = %f\n", seq_time);
	test();
  
  
  free(a);
}

/** -------------- SUB-PROCEDURES  ----------------- **/ 

/** procedure test() : verify sort results **/
void test() {
  int pass = 1;
  int i;
  for (i = 1; i < N; i++) {
    pass &= (a[i-1] <= a[i]);
  }

  printf(" TEST %s\n",(pass) ? "PASSed" : "FAILed");
}


/** procedure init() : initialize array "a" with data **/
void init() {
  int i;
  for (i = 0; i < N; i++) {
    a[i] = rand() % N; // (N - i);
  }
}

/** procedure  print() : print array elements **/
void print() {
  int i;
  for (i = 0; i < N; i++) {
    printf("%d\n", a[i]);
  }
  printf("\n");
}


/** INLINE procedure exchange() : pair swap **/
inline void exchange(int i, int j) {
  int t;
  t = a[i];
  a[i] = a[j];
  a[j] = t;
}



/** procedure compare() 
   The parameter dir indicates the sorting direction, ASCENDING 
   or DESCENDING; if (a[i] > a[j]) agrees with the direction, 
   then a[i] and a[j] are interchanged.
**/
inline void compare(int i, int j, int dir) {
  if (dir==(a[i]>a[j])) 
    exchange(i,j);
}

int CompAsc(const void *x,const void *y){
	const int *ix=(const int*)x;
	const int *iy=(const int*)y;
	return *ix-*iy;
}

int CompDesc(const void *x,const void *y){
	const int *ix=(const int*)x;
	const int *iy=(const int*)y;
	return *(iy)-*(ix);
	
}

/*------------------------------------------------------------------------------------------------------*/
/*******************************************************************************************************/
/*************PARALLEL CODE****************************************************************************/
/*****************************************************************************************************/
/*--------------------------------------------------------------------------------------------------*/

/** Procedure p_bitonicMerge() 
   It recursively sorts a bitonic sequence in ascending order, 
   if dir = ASCENDING, and in descending order otherwise. 
   The sequence to be sorted starts at index position lo,
   the parameter cbt is the number of elements to be sorted. Do it parallel.
 **/
void *p_bitonicMerge(void *ar) {

	int k,lo,hig,dire;
	element *arr=(element *)ar;		
	lo=arr->low;

	if (arr->high>1) {						//Checking if we are in leaves(bottom of tree)
											//if no,checking if active threads(cout1) are less 
		if(cout1<q){						//than maximum of active threads = q
			pthread_t thread0,thread1;
			element arr0,arr1;				//if yes,continue parallel sorting
			int rc2,rc3,i;
			k=(arr->high)/2;
			for (i=arr->low; i<(arr->low)+k; i++){
				compare(i, i+k, arr->dir);
			}
			arr0.low=arr->low;
			arr0.high=k;
			arr0.dir=arr->dir;
	 
			arr1.low=arr->low+k;
			arr1.high=k;
			arr1.dir=arr->dir;
			
			pthread_mutex_lock(&mute1);
			cout1+=2;							//Locking code for counter increase 
			pthread_mutex_unlock(&mute1);
			
			rc2=pthread_create(&thread0,NULL,p_bitonicMerge,&arr0);
			if(rc2){													//creating two threads
			  printf("\n ERROR : return code from create is...%d \n",rc2);  //one for each half
			}
		   
			rc3=pthread_create(&thread1,NULL,p_bitonicMerge,&arr1);
			if(rc3){													
			  printf("\n ERROR : return code from create is...%d \n",rc3);
			}
			
			pthread_join(thread0,NULL);
			pthread_join(thread1,NULL);
			
			pthread_mutex_lock(&mute);
			cout1=cout1-2;					//Locking code for counter decrease
			pthread_mutex_unlock(&mute);
			
		}
		else{
			
		bitonicMerge(arr->low,arr->high,arr->dir);		//if counter bigger than maximum threads(=q)
			return 0;										//do it serial
		}
  }

}
void *p_recBitonicSort(void *arr) {
	element *array=(element *)arr;
	int lo=array->low;
	int hig=array->high;
	int dire=array->dir;
	int i,mid,rc0,rc1,k;
	//pthread_t  self;	
	if((array->high)>1){

	if(cout1<q){	
		 k=(array->high)/2;
		
		 pthread_t thread0,thread1;		//thread0 for first half,thread1 for second
		 element matrix0,matrix1;		//same for matrix0,matrix1 .
		  
		  matrix0.low=array->low; 
		  matrix0.high=(array->high)/2;		//Next input for recursion	
		  matrix0.dir=ASCENDING;
		  matrix0.thd=cout1;
		  matrix1.low=((array->low)+k);
		  matrix1.high=(array->high)/2;	
		  matrix1.dir=DESCENDING;
		  matrix1.thd=cout1;

			pthread_mutex_lock(&mute);     //Lock the counter for active threads,increasing 
			cout1+=2;						//+2 for each MASTER thread
			pthread_mutex_unlock(&mute);
		  rc0=pthread_create(&thread0,NULL,p_recBitonicSort,&matrix0);
		  if(rc0){
			printf("\n ERROR : return code from create is...%d \n",rc0);
		  }
		  rc1=pthread_create(&thread1,NULL,p_recBitonicSort,&matrix1);
		  if(rc1){
			printf("\n ERROR : return code from create is...%d \n",rc1);
		  }
		  
		  pthread_join(thread0,NULL);
		  pthread_join(thread1,NULL);
		  pthread_mutex_lock(&mute);
		  cout1=cout1-2;				//After returning,decreace the counter 
		  pthread_mutex_unlock(&mute);  //we want the active number
		  //printf("calling merge from sort\n");
		 // array->high=(array->high)/2;
		element array1;
		array1.high=array->high;
		array1.low=array->low;
		array1.thd=cout1;				//Struct as input for merge
		array1.dir=array->dir;
		  p_bitonicMerge(&array1);		//Merging
	}
		   
   else{
	   
	 printf("enter case \n");	
//If counter is in upper bound=q(max threads)
	 
	recBitonicSort(array->low,array->high,array->dir); 	//then,do it serial

	}
 }
		
}

void p_sort() {
		element arr;
  
		arr.low=0;
		arr.high=N;
		arr.dir=ASCENDING;
		arr.thd=0;
		//printf("%d\n",arr.high);  //check.
		p_recBitonicSort(&arr);			//Start parallel bitonic sort ,passing struct arr 
		printf("sort over\n");
}

/***********************************SERIAL*********************************************************************/
/**************************************************************************************************************/



/** function recBitonicSort() 
    first produces a bitonic sequence by recursively sorting 
    its two halves in opposite sorting orders, and then
    calls bitonicMerge to make them in the same order
 **/
void bitonicMerge(int lo, int cnt, int dir) {
	
  if (cnt>1) {
	  
    int k=cnt/2;
    int i;
    for (i=lo; i<lo+k; i++)
      compare(i, i+k, dir);
    bitonicMerge(lo, k, dir);
    bitonicMerge(lo+k, k, dir);

  }
}


void recBitonicSort(int lo, int cnt, int dir) {
  if (cnt>1) {
    int k=cnt/2;
	//printf("%d\n",k);
    recBitonicSort(lo, k, ASCENDING);
    recBitonicSort(lo+k, k, DESCENDING);
    bitonicMerge(lo, cnt, dir);
  }
}


/** function sort() 
   Caller of recBitonicSort for sorting the entire array of length N 
   in ASCENDING order
 **/
void sort() {
  recBitonicSort(0, N, ASCENDING);
}



/*
  imperative version of bitonic sort
*/
void impBitonicSort() {

  int i,j,k;
  
  for (k=2; k<=N; k=2*k) {
    for (j=k>>1; j>0; j=j>>1) {
      for (i=0; i<N; i++) {
	int ij=i^j;
	if ((ij)>i) {
	  if ((i&k)==0 && a[i] > a[ij]) 
	      exchange(i,ij);
	  if ((i&k)!=0 && a[i] < a[ij])
	      exchange(i,ij);
	}
      }
    }
  }

}
