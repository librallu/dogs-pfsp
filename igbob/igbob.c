// main.c
// The file is self-contained, i.e. all functions required are described here. Only standard I/O and time libraries from C are required
// For doubts and questions, please send an email to:
// Victor Fernandez-Viagas - vfernandezviagas@us.es
// Jose M Framinan - framinan@us.es

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/timeb.h>

// define maximum number
#define MAX_LONG 2147483647

// define flowshop structure
typedef struct FLOWSHOP {

    int jobs;
    int machines;
    int ** pt;

} FLOWSHOP;

// declaration of functions employed in the file:
void brqsortLVectorD(long int *vector, int *index, int left, int right);
void calculate_e_q_f(int m,int k, int ** p_ij,int * partial_sequence,int next_job,long int ** e,long int ** q, long int ** f);
int Construction_FF(FLOWSHOP my_flowshop,int initial_length,int* partial_sequence,int * destructed_jobs);
void copyIVector(int *source, int *destination, int len);
void ** dyn_mat(int rows,int cols,int size);
int extractIVector(int * vector, int len, int pos);
void free_mat_int(int **m, int rows);
void free_mat_long(long int **m, int rows);
int IG_FF(FLOWSHOP my_flowshop,int d,double T,float stopping_criterion);
int IG_FF_outIter(FLOWSHOP my_flowshop,int d,double T,float stopping_criterion,int* iteraciones);
void insertIVector(int *vector, int length, int value, int pos);
long int lMax(long int a, long int b);
int ** loadPTimes_nrows(char *filename, int *jobs, int *machs);
int LS_FF(FLOWSHOP my_flowshop, int * initial_sequence,int length,int initial_makespan, float initial_time, float stopping_criterion, int * sequence_out);
int NEH_FF(FLOWSHOP my_flowshop);
int NEH_FF_seq(FLOWSHOP my_flowshop,int* partial_sequence);
int Pc_FF(FLOWSHOP my_flowshop, int * initial_sequence,int length,int initial_makespan, float initial_time, float stopping_criterion, int * sequence_out);
int LS(FLOWSHOP my_flowshop, int * initial_sequence,int length,int initial_makespan, float initial_time, float stopping_criterion, int * sequence_out);
void setval_Dvector(double * vector, int len, double val);
void setval_Ivector(int * vector, int len, int val);
void setval_Lmatrix(long int ** matrix, int rows, int cols, long int val);
void setval_Lvector(long * vector, int len, long int val);
int LS_FF(FLOWSHOP my_flowshop, int * initial_sequence,int length,int initial_makespan, float initial_time, float stopping_criterion, int * sequence_out);
int IterativeImprovement_Insertion(FLOWSHOP my_flowshop, int* initial_sequence,int length,int initial_makespan, float initial_time, float stopping_criterion, int* sequence_out);
int iteratedGreedy_Ruiz2007_Construction(FLOWSHOP my_flowshop,int initial_length,int* partial_sequence,int* destructed_jobs);

// functions:

void brqsortLVectorD(long int *vector, int *index, int left, int right) {


	register int i,j;
	long int x,y;
	int y_ind;

	i = left;
	j = right;
	x = vector[(left+right)/2];

	do {

		while( (vector[i]>x) && (i<right) ) i++;
		while( (x>vector[j]) && (j>left) ) j--;

		if(i<=j) {
			y = vector[i];
			y_ind = index[i];
			vector[i] = vector[j];
			index[i] = index[j];
			vector[j] = y;
			index[j] = y_ind;
			i++;
			j--;
		}

	} while(i<=j);

	if(left<j) brqsortLVectorD(vector, index, left, j);
	if(i<right) brqsortLVectorD(vector, index, i, right);

}   // end of brqsortLVectorD()

void copyIVector(int *source, int *destination, int len) {

	register int i;
	for(i=0;i<len;i++) *(destination+i)=*(source+i);

}	// end of copyIVector()


long int lMax(long int a, long int b) {

	if(a>=b) return a;
	return b;

}	// end of lMax()

void ** dyn_mat(int rows,int cols,int size)
 {

    register int i;
    void **m;

    // pointer to pointers' array
    m=(void**) malloc((unsigned)rows*sizeof(void*));

    for(i=0;i<rows;i++)
	{
	m[i]=malloc((unsigned)cols*size);
    }

return m;

}
	// end of dyn_mat()

void setval_Ivector(int * vector, int len, int val) {

	register int i;
	for(i=0;i<len;i++) vector[i] = val;

}	// end of setval_Ivector()

void setval_Lvector(long * vector, int len, long int val) {

	register int i;
	for(i=0;i<len;i++) vector[i] = val;

}	// end of setval_Lvector()

void setval_Lmatrix(long int ** matrix, int rows, int cols, long int val) {
	register int i,j;
	for(i=0;i<rows;i++) {
		for(j=0;j<cols;j++) matrix[i][j] = val;
	}

}   // end of setval_Lmatrix()

int extractIVector(int * vector, int len, int pos) {

    register int i;
    int * copy;

    // checking position is smaller than len

    copy = (int *) malloc(sizeof(int)*len);
    copyIVector(vector, copy, len);
    // displacement of all components from position pos
    for(i=pos;i< (len-1);i++) vector[i] = copy[i+1];

    free(copy);

    return (len-1);

}	// end of extractIVector()

void insertIVector(int *vector, int length, int value, int pos) {

	register int i;

	// shift forward those position downstream pos
	for(i=(length-1);i>pos;i--) *(vector+i)=*(vector+i-1);

	// set the value in position pos
	*(vector+pos)=value;

}		// end of insertIVector()

void free_mat_long(long int **m, int rows) {

	register int i;
	for(i=0;i<rows;i++) free(m[i]);
	free(m);

}

void free_mat_int(int **m, int rows) {

	register int i;
	for(i=0;i<rows;i++) free(m[i]);
	free(m);

}

 int ** loadPTimes_nrows(char *filename, int *jobs, int *machs) {

	// 0. Variables
	FILE *data;			// file pointer
	int temp_data;			// to store processing times
	register int i,j;		// counters
	 int ** pt;			// processing times matrix

	// 1. Checks that file exists
	if (!(data = fopen(filename,"rt"))) {
	    printf("File not found: ");
        exit(0);
    }

	// 2. Reading data...
	// 2.1. Headers
	fscanf(data,"%d\n", &temp_data);
	*(jobs) = temp_data;
	fscanf(data,"%d\n", &temp_data);
	*(machs) = temp_data;

	// 2.2. Ask for space for matrix pt
    pt = (int **) dyn_mat((*machs),(*jobs), sizeof(int));
	// 2.3. Processing times
	for(i=0;i<(*machs);i++) {
		for(j=0;j<*(jobs);j++) {
			fscanf(data,"%10d", &temp_data );
			// storing processing times
			pt[i][j] = temp_data;
		}
		fscanf(data,"\n");
	}	// end of 2.2.

	// 3. Close file
	fclose(data);

	// 5. Returning processing times matrix
        return pt;

}	// end of loadPTimes_nrows()

void setval_Dvector(double * vector, int len, double val) {

	register int i;
	for(i=0;i<len;i++) vector[i] = val;

}	// end of setval_Dvector()

void calculate_e_q_f(int m,int k, int ** p_ij,int * partial_sequence,int next_job,long int ** e,long int ** q, long int ** f)
{
    int i,j,l;
    // compute earliest starting time
    // first job, first machine
    e[0][0] = (long) p_ij[0][partial_sequence[0]];
    // first job, rest of the machines
    for(i=1;i< m; i++) {
        e[i][0] = e[i-1][0] + (long) p_ij[i][partial_sequence[0]];
    }
    // rest of the jobs
    for(j=1;j< k;j++) {
        // rest of the jobs, first machine
        e[0][j] = e[0][j-1] + (long) p_ij[0][partial_sequence[j]];
        // rest of the jobs, rest of the machines
        for(i=1;i< m;i++) {
            e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) p_ij[i][partial_sequence[j]];
        }
    }

    // compute tail
    // last job, last machine
    q[m-1][k-1] = (long) p_ij[m-1][partial_sequence[k-1]];
    // last job, rest of the machines
    for(i= m-2;i>=0; i--) {
        q[i][k-1] = q[i+1][k-1] + (long) p_ij[i][partial_sequence[k-1]];
    }
    // rest of the jobs
    for(j=k-2;j>=0;j--) {
        // rest of the jobs, last machine
        q[m-1][j] = q[m-1][j+1] + (long) p_ij[m-1][partial_sequence[j]];
        // rest of the jobs, rest of the machines
        for(i=m-2;i>=0;i--) {
            q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) p_ij[i][partial_sequence[j]];
        }
    }

    // compute f (for all positions where the job can be inserted)
    // first position (l=0)
    // first position, first machine
    f[0][0] = (long) p_ij[0][next_job];
    // first position, rest of the machines
    for(i=1;i< m;i++) {
        f[i][0] = f[i-1][0] + (long) p_ij[i][next_job];
    }
    // rest of the positions
    for(l=1;l<=k;l++) {
        // first machine
        f[0][l] = e[0][l-1] + (long) p_ij[0][next_job];
        // rest of the machines
        for(i=1;i<m;i++) {
            f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) p_ij[i][next_job];
        }
    }
}

int IG_FF(FLOWSHOP my_flowshop,int d,double T,float stopping_criterion)
{
	struct timeb init_time;
	struct timeb actual_time;
	float secs;
    ftime(&init_time);
    int * sequence_ini=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_iteration=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_best_local=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_destr_constr=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int mksp_initial=NEH_FF_seq(my_flowshop,sequence_ini);
    int i,j,l;
    int random_number;double random_double,temperature_aux;
    int stay_in_loop;int number_destructed_jobs;int is_destructed;
    long int best_makespan,makespan_destr_constr,best_local_makespan;
    long int makespan_iteration = MAX_LONG;
    int N=my_flowshop.jobs;
    int M=my_flowshop.machines;
     int ** p_ij = (int **) dyn_mat(M,N,sizeof(int));
    int * best_sequence=(int *) malloc(sizeof(int)*N);
    int * destruction_jobs=(int *) malloc(sizeof(int)*N-d);
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            p_ij[i][j]=my_flowshop.pt[i][j];
        }
    }
    double suma_p_ij = 0;
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < M; i++)
        {
            suma_p_ij = suma_p_ij +(double) p_ij[i][j];
        }
    }
    double Temperature = T * suma_p_ij / (double)N / (double)M / (double)10;
    ftime(&actual_time);
    secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);

    //LOCAL SEARCH
    best_local_makespan=LS_FF(my_flowshop,sequence_ini,N,mksp_initial,secs,stopping_criterion,sequence_iteration);
    best_makespan=best_local_makespan;
    copyIVector(sequence_iteration,best_sequence,N);
    copyIVector(sequence_iteration,sequence_best_local,N);
    //for(i=0;i<2;i++)
    while(secs<stopping_criterion)
    {
        //DESTRUCTION PHASE
        for (l = 0; l < d; l++)
        {
            destruction_jobs[l] = -1;
        }
        //Chosen random jobs
        for (l = 0; l < d; l++)
        {
            stay_in_loop = 1;
            random_number = rand()%N;
            while (stay_in_loop)
            {
                stay_in_loop = 0;
                random_number = rand()%N;
                for (i = 0; i < d && stay_in_loop == 0; i++)
                {
                    if (random_number == destruction_jobs[i])
                    {
                        stay_in_loop = 1;
                    }
                }
            }
            destruction_jobs[l] = random_number;
        }

        //CONSTRUCTION PHASE
        number_destructed_jobs = 0;
        for (j = 0; j < N; j++)
        {
            is_destructed =0;
            for (l = 0; l < d; l++)
            {
                if (sequence_iteration[j] == destruction_jobs[l])
                {
                    is_destructed = 1;
                }
            }
            if (is_destructed == 0)
            {
                sequence_destr_constr[j - number_destructed_jobs] = sequence_iteration[j];
            }
            else
            {
                number_destructed_jobs++;
            }
        }
        makespan_destr_constr=Construction_FF(my_flowshop,N-d,sequence_destr_constr,destruction_jobs);

        //LOCAL SEARCH
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
        best_local_makespan=LS_FF(my_flowshop,sequence_destr_constr,N,makespan_destr_constr,secs,stopping_criterion,sequence_best_local);
        random_double=(double) rand()/(double) RAND_MAX;
        if(best_local_makespan<makespan_iteration)
        {
            makespan_iteration=best_local_makespan;
            copyIVector(sequence_best_local,sequence_iteration,N);
            if(best_local_makespan<best_makespan)
            {
                best_makespan=best_local_makespan;
                copyIVector(sequence_best_local,best_sequence,N);
            }
        }
        else
        {
            temperature_aux=-(double)(best_local_makespan - makespan_iteration) / Temperature;
            if (random_double <=  exp(temperature_aux))
            {
                makespan_iteration = best_local_makespan;
                copyIVector(sequence_best_local,sequence_iteration,N);
            }
        }
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    }
    free(sequence_best_local);
    free(sequence_ini);
    free(sequence_destr_constr);
    free(destruction_jobs);
    free(sequence_iteration);
    free(best_sequence);
    free_mat_int(p_ij,M);
    return best_makespan;
}

int IG_FF_outIter(FLOWSHOP my_flowshop,int d,double T,float stopping_criterion,int* iteraciones)
{
	struct timeb init_time;
	struct timeb actual_time;
	float secs;
    ftime(&init_time);
    int * sequence_ini=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_iteration=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_best_local=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_destr_constr=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int mksp_initial=NEH_FF_seq(my_flowshop,sequence_ini);
    int i,j,l;
    int random_number;double random_double,temperature_aux;
    int stay_in_loop;int number_destructed_jobs;int is_destructed;
    long int best_makespan, makespan_destr_constr,best_local_makespan;
    long int makespan_iteration = MAX_LONG;
    int N=my_flowshop.jobs;
    int M=my_flowshop.machines;
     int ** p_ij = (int **) dyn_mat(M,N,sizeof(int));
    int * best_sequence=(int *) malloc(sizeof(int)*N);
    int * destruction_jobs=(int *) malloc(sizeof(int)*N-d);
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            p_ij[i][j]=my_flowshop.pt[i][j];
        }
    }
    double suma_p_ij = 0;
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < M; i++)
        {
            suma_p_ij = suma_p_ij +(double) p_ij[i][j];
        }
    }
    double Temperature = T * suma_p_ij / (double)N / (double)M / (double)10;
    ftime(&actual_time);
    secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);


    //LOCAL SEARCH
    best_local_makespan=LS_FF(my_flowshop,sequence_ini,N,mksp_initial,secs,stopping_criterion,sequence_iteration);
    best_makespan=best_local_makespan;
    copyIVector(sequence_iteration,best_sequence,N);
    copyIVector(sequence_iteration,sequence_best_local,N);
    *iteraciones=0;
    while(secs<stopping_criterion)
    {
        *iteraciones=*iteraciones+1;
        //DESTRUCTION PHASE
        for (l = 0; l < d; l++)
        {
            destruction_jobs[l] = -1;
        }
        //Chosen random jobs
        for (l = 0; l < d; l++)
        {
            stay_in_loop = 1;
            random_number = rand()%N;
            while (stay_in_loop)
            {
                stay_in_loop = 0;
                random_number = rand()%N;
                for (i = 0; i < d && stay_in_loop == 0; i++)
                {
                    if (random_number == destruction_jobs[i])
                    {
                        stay_in_loop = 1;
                    }
                }
            }
            destruction_jobs[l] = random_number;
        }

        //CONSTRUCTION PHASE
        number_destructed_jobs = 0;
        for (j = 0; j < N; j++)
        {
            is_destructed =0;
            for (l = 0; l < d; l++)
            {
                if (sequence_iteration[j] == destruction_jobs[l])
                {
                    is_destructed = 1;
                }
            }
            if (is_destructed == 0)
            {
                sequence_destr_constr[j - number_destructed_jobs] = sequence_iteration[j];
            }
            else
            {
                number_destructed_jobs++;
            }
        }
        makespan_destr_constr=Construction_FF(my_flowshop,N-d,sequence_destr_constr,destruction_jobs);

        //LOCAL SEARCH
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
        best_local_makespan=LS_FF(my_flowshop,sequence_destr_constr,N,makespan_destr_constr,secs,stopping_criterion,sequence_best_local);
        random_double=(double) rand()/(double) RAND_MAX;
        if(best_local_makespan<makespan_iteration)
        {
            makespan_iteration=best_local_makespan;
            copyIVector(sequence_best_local,sequence_iteration,N);
            if(best_local_makespan<best_makespan)
            {
                best_makespan=best_local_makespan;
                copyIVector(sequence_best_local,best_sequence,N);
            }
        }
        else
        {
            temperature_aux=-(double)(best_local_makespan - makespan_iteration) / Temperature;
            if (random_double <= exp(temperature_aux))
            {
                makespan_iteration = best_local_makespan;
                copyIVector(sequence_best_local,sequence_iteration,N);
            }
        }
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    }
    free(sequence_best_local);
    free(sequence_ini);
    free(sequence_destr_constr);
    free(destruction_jobs);
    free(sequence_iteration);
    free(best_sequence);
    free_mat_int(p_ij,M);
    return best_makespan;
}


int IG_outIter(FLOWSHOP my_flowshop,int d,double T,float stopping_criterion,int* iteraciones)
{
	struct timeb init_time;
	struct timeb actual_time;
	float secs;
    ftime(&init_time);
    int * sequence_ini=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_iteration=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_best_local=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_destr_constr=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int mksp_initial=NEH_seq(my_flowshop,sequence_ini);
    int i,j,l;
    int random_number;double random_double,temperature_aux;
    int stay_in_loop;int number_destructed_jobs;int is_destructed;
    long int best_makespan, makespan_destr_constr,best_local_makespan;
    long int makespan_iteration = MAX_LONG;
    int N=my_flowshop.jobs;
    int M=my_flowshop.machines;
     int ** p_ij = (int **) dyn_mat(M,N,sizeof(int));
    int * best_sequence=(int *) malloc(sizeof(int)*N);
    int * destruction_jobs=(int *) malloc(sizeof(int)*N-d);
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            p_ij[i][j]=my_flowshop.pt[i][j];
        }
    }
    double suma_p_ij = 0;
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < M; i++)
        {
            suma_p_ij = suma_p_ij +(double) p_ij[i][j];
        }
    }
    double Temperature = T * suma_p_ij / (double)N / (double)M / (double)10;
    ftime(&actual_time);
    secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);


    //LOCAL SEARCH
    best_local_makespan=LS(my_flowshop,sequence_ini,N,mksp_initial,secs,stopping_criterion,sequence_iteration);
    best_makespan=best_local_makespan;
    copyIVector(sequence_iteration,best_sequence,N);
    copyIVector(sequence_iteration,sequence_best_local,N);
    *iteraciones=0;
    while(secs<stopping_criterion)
    {
        *iteraciones=*iteraciones+1;
        //DESTRUCTION PHASE
        for (l = 0; l < d; l++)
        {
            destruction_jobs[l] = -1;
        }
        //Chosen random jobs
        for (l = 0; l < d; l++)
        {
            stay_in_loop = 1;
            random_number = rand()%N;
            while (stay_in_loop)
            {
                stay_in_loop = 0;
                random_number = rand()%N;
                for (i = 0; i < d && stay_in_loop == 0; i++)
                {
                    if (random_number == destruction_jobs[i])
                    {
                        stay_in_loop = 1;
                    }
                }
            }
            destruction_jobs[l] = random_number;
        }

        //CONSTRUCTION PHASE
        number_destructed_jobs = 0;
        for (j = 0; j < N; j++)
        {
            is_destructed =0;
            for (l = 0; l < d; l++)
            {
                if (sequence_iteration[j] == destruction_jobs[l])
                {
                    is_destructed = 1;
                }
            }
            if (is_destructed == 0)
            {
                sequence_destr_constr[j - number_destructed_jobs] = sequence_iteration[j];
            }
            else
            {
                number_destructed_jobs++;
            }
        }
        makespan_destr_constr=Construction(my_flowshop,N-d,sequence_destr_constr,destruction_jobs);

        //LOCAL SEARCH
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
        best_local_makespan=LS(my_flowshop,sequence_destr_constr,N,makespan_destr_constr,secs,stopping_criterion,sequence_best_local);
        random_double=(double) rand()/(double) RAND_MAX;
        if(best_local_makespan<makespan_iteration)
        {
            makespan_iteration=best_local_makespan;
            copyIVector(sequence_best_local,sequence_iteration,N);
            if(best_local_makespan<best_makespan)
            {
                best_makespan=best_local_makespan;
                copyIVector(sequence_best_local,best_sequence,N);
            }
        }
        else
        {
            temperature_aux=-(double)(best_local_makespan - makespan_iteration) / Temperature;
            if (random_double <= exp(temperature_aux))
            {
                makespan_iteration = best_local_makespan;
                copyIVector(sequence_best_local,sequence_iteration,N);
            }
        }
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    }
    free(sequence_best_local);
    free(sequence_ini);
    free(sequence_destr_constr);
    free(destruction_jobs);
    free(sequence_iteration);
    free(best_sequence);
    free_mat_int(p_ij,M);
    return best_makespan;
}

int NEH_FF_seq(FLOWSHOP my_flowshop,int* partial_sequence) {

    register int i,j,k,l;
    long int best_of_position_l;
    long int curr_best = MAX_LONG;
    long int it_bp=MAX_LONG;
    long int idle_time;

    int bp = -1;
    int tb=0;
    int index=0;

    long * sum_processing_times = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * non_scheduled_jobs = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    long * array_makespan = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * ptb = (int *) malloc(sizeof(int)*my_flowshop.jobs);

    setval_Lvector(sum_processing_times,my_flowshop.jobs, 0);

    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);

    long * f_prima = (long *) malloc(sizeof(long)*my_flowshop.machines);
    setval_Lvector(f_prima,my_flowshop.machines,0);


    // store the sum of the completion times
    for(j=0;j<my_flowshop.jobs; j++) {
        for(i=0;i< my_flowshop.machines;i++) {
            sum_processing_times[j] += (long) my_flowshop.pt[i][j];
        }
    }

    // sort the jobs in descending order of the sum of the processing times
    for(j=0;j< my_flowshop.jobs;j++) non_scheduled_jobs[j] = j;
    brqsortLVectorD(sum_processing_times, non_scheduled_jobs, 0, my_flowshop.jobs-1);


    // start with the initial job (to be inserted in pos 0) and remove the job
    int tamano_actual=my_flowshop.jobs;
    partial_sequence[0] = non_scheduled_jobs[0];
    extractIVector(non_scheduled_jobs,my_flowshop.jobs,0);
    tamano_actual--;

    // for the rest of the jobs
    for(k=1;k < my_flowshop.jobs; k++) {

            // compute earliest starting time
            // first job, first machine
            e[0][0] = (long) my_flowshop.pt[0][partial_sequence[0]];
            // first job, rest of the machines
            for(i=1;i< my_flowshop.machines; i++) {
                e[i][0] = e[i-1][0] + (long) my_flowshop.pt[i][partial_sequence[0]];
            }
            // rest of the jobs
            for(j=1;j< k;j++) {
                // rest of the jobs, first machine
                e[0][j] = e[0][j-1] + (long) my_flowshop.pt[0][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=1;i< my_flowshop.machines;i++) {
                    e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute tail
            // last job, last machine
            q[my_flowshop.machines-1][k-1] = (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[k-1]];
            // last job, rest of the machines
            for(i= my_flowshop.machines-2;i>=0; i--) {
                q[i][k-1] = q[i+1][k-1] + (long) my_flowshop.pt[i][partial_sequence[k-1]];
            }
            // rest of the jobs
            for(j=k-2;j>=0;j--) {
                // rest of the jobs, last machine
                q[my_flowshop.machines-1][j] = q[my_flowshop.machines-1][j+1] + (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=my_flowshop.machines-2;i>=0;i--) {
                    q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute f (for all positions where the job can be inserted)
            // first position (l=0)
            // first position, first machine
            f[0][0] = (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
            // first position, rest of the machines
            for(i=1;i< my_flowshop.machines;i++) {
                f[i][0] = f[i-1][0] + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
            }
            // rest of the positions
            for(l=1;l<=k;l++) {
                // first machine
                f[0][l] = e[0][l-1] + (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
                // rest of the machines
                for(i=1;i<my_flowshop.machines;i++) {
                    f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
                }
            }

            // compute best cmax
            curr_best = MAX_LONG;
            for(l=0;l<=k;l++)
            {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++)
                {
                    if( f[i][l] + q[i][l] > best_of_position_l)
                    best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < curr_best )
                {
                    curr_best = best_of_position_l;
                    bp = l;
                }
            }

            //number of tiesint contador = 0;
            tb=0;
            for(l=0;l<=k;l++)
            {
                if (array_makespan[l] == curr_best)
                {
                    tb++;
                }
            }

            //More than 1 tie
            if(tb>1&& k<my_flowshop.jobs-1)
            {
                index=0;
                //Find the position with ties
                for(l=0;l<=k;l++)
                {
                    if(array_makespan[l] == curr_best)
                    {
                        ptb[index]=l;
                        index++;
                    }
                }
                it_bp = MAX_LONG;
                for(j=0;j<tb;j++)
                {
                    idle_time=0;
                    if(ptb[j]<k)
                    {
                        f_prima[0]=f[0][ptb[j]]+(long)my_flowshop.pt[0][partial_sequence[ptb[j]]];
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            f_prima[i]=lMax(f_prima[i-1],f[i][ptb[j]])+(long)my_flowshop.pt[i][partial_sequence[ptb[j]]];
                            idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]]+(long)my_flowshop.pt[i][partial_sequence[ptb[j]]]+lMax(0,f_prima[i-1]-f[i][ptb[j]]);
                        }
                    }
                    else
                    {
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]-1];
                        }
                    }
                    if(idle_time<it_bp)
                    {
                        it_bp=idle_time;
                        bp=ptb[j];
                    }
                }
            }

            // inserts the job in the best position
            insertIVector(partial_sequence,k+1,non_scheduled_jobs[0],bp);

            // removes the job from the non_scheduled
            extractIVector(non_scheduled_jobs,tamano_actual,0);
            tamano_actual--;

    }

    free(non_scheduled_jobs);
    free(sum_processing_times);
    free(array_makespan);
    free(ptb);
    free(f_prima);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);

    return curr_best;
}


int NEH_seq(FLOWSHOP my_flowshop,int* partial_sequence) {

    register int i,j,k,l;
    long int best_of_position_l;
    long int curr_best = MAX_LONG;
    long int it_bp=MAX_LONG;
    long int idle_time;

    int bp = -1;
    int tb=0;
    int index=0;

    long * sum_processing_times = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * non_scheduled_jobs = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    long * array_makespan = (long *) malloc(sizeof(long)*my_flowshop.jobs);

    setval_Lvector(sum_processing_times,my_flowshop.jobs, 0);

    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);

    // store the sum of the completion times
    for(j=0;j<my_flowshop.jobs; j++) {
        for(i=0;i< my_flowshop.machines;i++) {
            sum_processing_times[j] += (long) my_flowshop.pt[i][j];
        }
    }

    // sort the jobs in descending order of the sum of the processing times
    for(j=0;j< my_flowshop.jobs;j++) non_scheduled_jobs[j] = j;
    brqsortLVectorD(sum_processing_times, non_scheduled_jobs, 0, my_flowshop.jobs-1);


    // start with the initial job (to be inserted in pos 0) and remove the job
    int tamano_actual=my_flowshop.jobs;
    partial_sequence[0] = non_scheduled_jobs[0];
    extractIVector(non_scheduled_jobs,my_flowshop.jobs,0);
    tamano_actual--;

    // for the rest of the jobs
    for(k=1;k < my_flowshop.jobs; k++) {

            // compute earliest starting time
            // first job, first machine
            e[0][0] = (long) my_flowshop.pt[0][partial_sequence[0]];
            // first job, rest of the machines
            for(i=1;i< my_flowshop.machines; i++) {
                e[i][0] = e[i-1][0] + (long) my_flowshop.pt[i][partial_sequence[0]];
            }
            // rest of the jobs
            for(j=1;j< k;j++) {
                // rest of the jobs, first machine
                e[0][j] = e[0][j-1] + (long) my_flowshop.pt[0][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=1;i< my_flowshop.machines;i++) {
                    e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute tail
            // last job, last machine
            q[my_flowshop.machines-1][k-1] = (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[k-1]];
            // last job, rest of the machines
            for(i= my_flowshop.machines-2;i>=0; i--) {
                q[i][k-1] = q[i+1][k-1] + (long) my_flowshop.pt[i][partial_sequence[k-1]];
            }
            // rest of the jobs
            for(j=k-2;j>=0;j--) {
                // rest of the jobs, last machine
                q[my_flowshop.machines-1][j] = q[my_flowshop.machines-1][j+1] + (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=my_flowshop.machines-2;i>=0;i--) {
                    q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute f (for all positions where the job can be inserted)
            // first position (l=0)
            // first position, first machine
            f[0][0] = (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
            // first position, rest of the machines
            for(i=1;i< my_flowshop.machines;i++) {
                f[i][0] = f[i-1][0] + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
            }
            // rest of the positions
            for(l=1;l<=k;l++) {
                // first machine
                f[0][l] = e[0][l-1] + (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
                // rest of the machines
                for(i=1;i<my_flowshop.machines;i++) {
                    f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
                }
            }

            // compute best cmax
            curr_best = MAX_LONG;
            for(l=0;l<=k;l++)
            {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++)
                {
                    if( f[i][l] + q[i][l] > best_of_position_l)
                    best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < curr_best )
                {
                    curr_best = best_of_position_l;
                    bp = l;
                }
            }

            // inserts the job in the best position
            insertIVector(partial_sequence,k+1,non_scheduled_jobs[0],bp);

            // removes the job from the non_scheduled
            extractIVector(non_scheduled_jobs,tamano_actual,0);
            tamano_actual--;

    }

    free(non_scheduled_jobs);
    free(sum_processing_times);
    free(array_makespan);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);

    return curr_best;
}


int NEH_FF(FLOWSHOP my_flowshop) {

    register int i,j,k,l;
    long int best_of_position_l;
    long int curr_best = MAX_LONG;
    long int it_bp=MAX_LONG;
    long int idle_time;

    int bp = -1;
    int tb=0;
    int index=0;

    int * partial_sequence = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    setval_Ivector(partial_sequence,my_flowshop.jobs, 0);
    long * sum_processing_times = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * non_scheduled_jobs = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    long * array_makespan = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * ptb = (int *) malloc(sizeof(int)*my_flowshop.jobs);

    setval_Lvector(sum_processing_times,my_flowshop.jobs, 0);

    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);

    long * f_prima = (long *) malloc(sizeof(long)*my_flowshop.machines);
    setval_Lvector(f_prima,my_flowshop.machines,0);


    // store the sum of the completion times
    for(j=0;j<my_flowshop.jobs; j++) {
        for(i=0;i< my_flowshop.machines;i++) {
            sum_processing_times[j] += (long) my_flowshop.pt[i][j];
        }
    }

    // sort the jobs in descending order of the sum of the processing times
    for(j=0;j< my_flowshop.jobs;j++) non_scheduled_jobs[j] = j;
    brqsortLVectorD(sum_processing_times, non_scheduled_jobs, 0, my_flowshop.jobs-1);

    // start with the initial job (to be inserted in pos 0) and remove the job
    int tamano_actual=my_flowshop.jobs;
    partial_sequence[0] = non_scheduled_jobs[0];
    extractIVector(non_scheduled_jobs,my_flowshop.jobs,0);
    tamano_actual--;

    // for the rest of the jobs
    for(k=1;k < my_flowshop.jobs; k++) {

            // compute earliest starting time
            // first job, first machine
            e[0][0] = (long) my_flowshop.pt[0][partial_sequence[0]];
            // first job, rest of the machines
            for(i=1;i< my_flowshop.machines; i++) {
                e[i][0] = e[i-1][0] + (long) my_flowshop.pt[i][partial_sequence[0]];
            }
            // rest of the jobs
            for(j=1;j< k;j++) {
                // rest of the jobs, first machine
                e[0][j] = e[0][j-1] + (long) my_flowshop.pt[0][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=1;i< my_flowshop.machines;i++) {
                    e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute tail
            // last job, last machine
            q[my_flowshop.machines-1][k-1] = (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[k-1]];
            // last job, rest of the machines
            for(i= my_flowshop.machines-2;i>=0; i--) {
                q[i][k-1] = q[i+1][k-1] + (long) my_flowshop.pt[i][partial_sequence[k-1]];
            }
            // rest of the jobs
            for(j=k-2;j>=0;j--) {
                // rest of the jobs, last machine
                q[my_flowshop.machines-1][j] = q[my_flowshop.machines-1][j+1] + (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=my_flowshop.machines-2;i>=0;i--) {
                    q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute f (for all positions where the job can be inserted)
            // first position (l=0)
            // first position, first machine
            f[0][0] = (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
            // first position, rest of the machines
            for(i=1;i< my_flowshop.machines;i++) {
                f[i][0] = f[i-1][0] + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
            }
            // rest of the positions
            for(l=1;l<=k;l++) {
                // first machine
                f[0][l] = e[0][l-1] + (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
                // rest of the machines
                for(i=1;i<my_flowshop.machines;i++) {
                    f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
                }
            }

            // compute best cmax
            curr_best = MAX_LONG;
            for(l=0;l<=k;l++)
            {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++)
                {
                    if( f[i][l] + q[i][l] > best_of_position_l)
                    best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < curr_best )
                {
                    curr_best = best_of_position_l;
                    bp = l;
                }
            }

            //number of tiesint contador = 0;
            tb=0;
            for(l=0;l<=k;l++)
            {
                if (array_makespan[l] == curr_best)
                {
                    tb++;
                }
            }

            //More than 1 tie
            if(tb>1&& k<my_flowshop.jobs-1)
            {
                index=0;
                //Find the position with ties
                for(l=0;l<=k;l++)
                {
                    if(array_makespan[l] == curr_best)
                    {
                        ptb[index]=l;
                        index++;
                    }
                }
                it_bp = MAX_LONG;
                for(j=0;j<tb;j++)
                {
                    idle_time=0;
                    if(ptb[j]<k)
                    {
                        f_prima[0]=f[0][ptb[j]]+(long)my_flowshop.pt[0][partial_sequence[ptb[j]]];
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            f_prima[i]=lMax(f_prima[i-1],f[i][ptb[j]])+(long)my_flowshop.pt[i][partial_sequence[ptb[j]]];
                            idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]]+(long)my_flowshop.pt[i][partial_sequence[ptb[j]]]+lMax(0,f_prima[i-1]-f[i][ptb[j]]);
                        }
                    }
                    else
                    {
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]-1];
                        }
                    }
                    if(idle_time<it_bp)
                    {
                        it_bp=idle_time;
                        bp=ptb[j];
                    }
                }
            }

            // inserts the job in the best position
            insertIVector(partial_sequence,k+1,non_scheduled_jobs[0],bp);

            // removes the job from the non_scheduled
            extractIVector(non_scheduled_jobs,tamano_actual,0);
            tamano_actual--;

    }

    free(non_scheduled_jobs);
    free(sum_processing_times);
    free(array_makespan);
    free(ptb);
    free(partial_sequence);
    free(f_prima);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);

    return curr_best;
}


int NEH(FLOWSHOP my_flowshop) {

    register int i,j,k,l;
    long int best_of_position_l;
    long int curr_best = MAX_LONG;
    long int it_bp=MAX_LONG;
    long int idle_time;

    int bp = -1;
    int tb=0;
    int index=0;

    int * partial_sequence = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    setval_Ivector(partial_sequence,my_flowshop.jobs, 0);
    long * sum_processing_times = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * non_scheduled_jobs = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    long * array_makespan = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    setval_Lvector(sum_processing_times,my_flowshop.jobs, 0);

    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);


    // store the sum of the completion times
    for(j=0;j<my_flowshop.jobs; j++) {
        for(i=0;i< my_flowshop.machines;i++) {
            sum_processing_times[j] += (long) my_flowshop.pt[i][j];
        }
    }

    // sort the jobs in descending order of the sum of the processing times
    for(j=0;j< my_flowshop.jobs;j++) non_scheduled_jobs[j] = j;
    brqsortLVectorD(sum_processing_times, non_scheduled_jobs, 0, my_flowshop.jobs-1);

    // start with the initial job (to be inserted in pos 0) and remove the job
    int tamano_actual=my_flowshop.jobs;
    partial_sequence[0] = non_scheduled_jobs[0];
    extractIVector(non_scheduled_jobs,my_flowshop.jobs,0);
    tamano_actual--;

    // for the rest of the jobs
    for(k=1;k < my_flowshop.jobs; k++) {

            // compute earliest starting time
            // first job, first machine
            e[0][0] = (long) my_flowshop.pt[0][partial_sequence[0]];
            // first job, rest of the machines
            for(i=1;i< my_flowshop.machines; i++) {
                e[i][0] = e[i-1][0] + (long) my_flowshop.pt[i][partial_sequence[0]];
            }
            // rest of the jobs
            for(j=1;j< k;j++) {
                // rest of the jobs, first machine
                e[0][j] = e[0][j-1] + (long) my_flowshop.pt[0][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=1;i< my_flowshop.machines;i++) {
                    e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute tail
            // last job, last machine
            q[my_flowshop.machines-1][k-1] = (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[k-1]];
            // last job, rest of the machines
            for(i= my_flowshop.machines-2;i>=0; i--) {
                q[i][k-1] = q[i+1][k-1] + (long) my_flowshop.pt[i][partial_sequence[k-1]];
            }
            // rest of the jobs
            for(j=k-2;j>=0;j--) {
                // rest of the jobs, last machine
                q[my_flowshop.machines-1][j] = q[my_flowshop.machines-1][j+1] + (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=my_flowshop.machines-2;i>=0;i--) {
                    q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute f (for all positions where the job can be inserted)
            // first position (l=0)
            // first position, first machine
            f[0][0] = (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
            // first position, rest of the machines
            for(i=1;i< my_flowshop.machines;i++) {
                f[i][0] = f[i-1][0] + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
            }
            // rest of the positions
            for(l=1;l<=k;l++) {
                // first machine
                f[0][l] = e[0][l-1] + (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
                // rest of the machines
                for(i=1;i<my_flowshop.machines;i++) {
                    f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
                }
            }

            // compute best cmax
            curr_best = MAX_LONG;
            for(l=0;l<=k;l++)
            {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++)
                {
                    if( f[i][l] + q[i][l] > best_of_position_l)
                    best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < curr_best )
                {
                    curr_best = best_of_position_l;
                    bp = l;
                }
            }

            // inserts the job in the best position
            insertIVector(partial_sequence,k+1,non_scheduled_jobs[0],bp);

            // removes the job from the non_scheduled
            extractIVector(non_scheduled_jobs,tamano_actual,0);
            tamano_actual--;

    }

    free(non_scheduled_jobs);
    free(sum_processing_times);
    free(array_makespan);
    free(partial_sequence);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);

    return curr_best;
}

int Construction_FF(FLOWSHOP my_flowshop,int initial_length,int* partial_sequence,int * destructed_jobs)
{

    register int i,j,k,l;
    long int best_of_position_l;
    long int curr_best = MAX_LONG;

    int bp = -1;

    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);

    long int it_bp=MAX_LONG;
    long int idle_time;
    int tb=0;
    int index=0;
    long * array_makespan = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * ptb = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    long * f_prima = (long *) malloc(sizeof(long)*my_flowshop.machines);
    setval_Lvector(f_prima,my_flowshop.machines,0);

    // for the rest of the jobs
    for(k=initial_length;k < my_flowshop.jobs; k++) {

            // compute earliest starting time
            // first job, first machine
            e[0][0] = (long) my_flowshop.pt[0][partial_sequence[0]];
            // first job, rest of the machines
            for(i=1;i< my_flowshop.machines; i++) {
                e[i][0] = e[i-1][0] + (long) my_flowshop.pt[i][partial_sequence[0]];
            }
            // rest of the jobs
            for(j=1;j< k;j++) {
                // rest of the jobs, first machine
                e[0][j] = e[0][j-1] + (long) my_flowshop.pt[0][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=1;i< my_flowshop.machines;i++) {
                    e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute tail
            // last job, last machine
            q[my_flowshop.machines-1][k-1] = (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[k-1]];
            // last job, rest of the machines
            for(i= my_flowshop.machines-2;i>=0; i--) {
                q[i][k-1] = q[i+1][k-1] + (long) my_flowshop.pt[i][partial_sequence[k-1]];
            }
            // rest of the jobs
            for(j=k-2;j>=0;j--) {
                // rest of the jobs, last machine
                q[my_flowshop.machines-1][j] = q[my_flowshop.machines-1][j+1] + (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=my_flowshop.machines-2;i>=0;i--) {
                    q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute f (for all positions where the job can be inserted)
            // first position (l=0)
            // first position, first machine
            f[0][0] = (long) my_flowshop.pt[0][destructed_jobs[k-initial_length]];
            // first position, rest of the machines
            for(i=1;i< my_flowshop.machines;i++) {
                f[i][0] = f[i-1][0] + (long) my_flowshop.pt[i][destructed_jobs[k-initial_length]];
            }
            // rest of the positions
            for(l=1;l<=k;l++) {
                // first machine
                f[0][l] = e[0][l-1] + (long) my_flowshop.pt[0][destructed_jobs[k-initial_length]];
                // rest of the machines
                for(i=1;i<my_flowshop.machines;i++) {
                    f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) my_flowshop.pt[i][destructed_jobs[k-initial_length]];
                }
            }

            // compute best cmax
            curr_best = MAX_LONG;
            for(l=0;l<=k;l++) {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++) {
                    if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < curr_best ) {
                    curr_best = best_of_position_l;
                    bp = l;
                }
            }


            //number of tiesint contador = 0;
            tb=0;
            for(l=0;l<=k;l++)
            {
                if (array_makespan[l] == curr_best)
                {
                    tb++;
                }
            }

            //More than 1 tie
            if(tb>1)
            {
                index=0;
                //Find the position with ties
                for(l=0;l<=k;l++)
                {
                    if(array_makespan[l] == curr_best)
                    {
                        ptb[index]=l;
                        index++;
                    }
                }
                it_bp = MAX_LONG;
                for(j=0;j<tb;j++)
                {
                    idle_time=0;
                    if(ptb[j]<k)
                    {
                        f_prima[0]=f[0][ptb[j]]+(long)my_flowshop.pt[0][partial_sequence[ptb[j]]];
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            f_prima[i]=lMax(f_prima[i-1],f[i][ptb[j]])+(long)my_flowshop.pt[i][partial_sequence[ptb[j]]];
                            idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]]+(long)my_flowshop.pt[i][partial_sequence[ptb[j]]]+lMax(0,f_prima[i-1]-f[i][ptb[j]]);
                        }
                    }
                    else
                    {
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]-1];
                        }
                    }
                    if(idle_time<it_bp)
                    {
                        it_bp=idle_time;
                        bp=ptb[j];
                    }
                }
            }

            // inserts the job in the best position
            insertIVector(partial_sequence,k+1,destructed_jobs[k-initial_length],bp);


    }


    free(ptb);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);
    free(array_makespan);
    free(f_prima);

    return curr_best;
}


int Construction(FLOWSHOP my_flowshop,int initial_length,int* partial_sequence,int * destructed_jobs)
{

    register int i,j,k,l;
    long int best_of_position_l;
    long int curr_best = MAX_LONG;

    int bp = -1;

    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);

    long int it_bp=MAX_LONG;
    long int idle_time;
    int tb=0;
    int index=0;
    long * array_makespan = (long *) malloc(sizeof(long)*my_flowshop.jobs);

    // for the rest of the jobs
    for(k=initial_length;k < my_flowshop.jobs; k++) {

            // compute earliest starting time
            // first job, first machine
            e[0][0] = (long) my_flowshop.pt[0][partial_sequence[0]];
            // first job, rest of the machines
            for(i=1;i< my_flowshop.machines; i++) {
                e[i][0] = e[i-1][0] + (long) my_flowshop.pt[i][partial_sequence[0]];
            }
            // rest of the jobs
            for(j=1;j< k;j++) {
                // rest of the jobs, first machine
                e[0][j] = e[0][j-1] + (long) my_flowshop.pt[0][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=1;i< my_flowshop.machines;i++) {
                    e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute tail
            // last job, last machine
            q[my_flowshop.machines-1][k-1] = (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[k-1]];
            // last job, rest of the machines
            for(i= my_flowshop.machines-2;i>=0; i--) {
                q[i][k-1] = q[i+1][k-1] + (long) my_flowshop.pt[i][partial_sequence[k-1]];
            }
            // rest of the jobs
            for(j=k-2;j>=0;j--) {
                // rest of the jobs, last machine
                q[my_flowshop.machines-1][j] = q[my_flowshop.machines-1][j+1] + (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=my_flowshop.machines-2;i>=0;i--) {
                    q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute f (for all positions where the job can be inserted)
            // first position (l=0)
            // first position, first machine
            f[0][0] = (long) my_flowshop.pt[0][destructed_jobs[k-initial_length]];
            // first position, rest of the machines
            for(i=1;i< my_flowshop.machines;i++) {
                f[i][0] = f[i-1][0] + (long) my_flowshop.pt[i][destructed_jobs[k-initial_length]];
            }
            // rest of the positions
            for(l=1;l<=k;l++) {
                // first machine
                f[0][l] = e[0][l-1] + (long) my_flowshop.pt[0][destructed_jobs[k-initial_length]];
                // rest of the machines
                for(i=1;i<my_flowshop.machines;i++) {
                    f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) my_flowshop.pt[i][destructed_jobs[k-initial_length]];
                }
            }

            // compute best cmax
            curr_best = MAX_LONG;
            for(l=0;l<=k;l++) {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++) {
                    if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < curr_best ) {
                    curr_best = best_of_position_l;
                    bp = l;
                }
            }

            // inserts the job in the best position
            insertIVector(partial_sequence,k+1,destructed_jobs[k-initial_length],bp);


    }


    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);
    free(array_makespan);

    return curr_best;
}



int LS(FLOWSHOP my_flowshop, int * initial_sequence,int length,int initial_makespan, float initial_time, float stopping_criterion, int * sequence_out)
{
    int i,j,k,l;
	struct timeb init_time;
	struct timeb actual_time;
	float secs;
    ftime(&init_time);
    int index,rnd_numb,temp,j_aux,job_chosen;
    int bp = -1;
    long int best_of_position_l;
    long int best_makespan = initial_makespan;
    long int min_local_mkspan=0;
    int * sub_seq = (int *) malloc(sizeof(int)*length);
    int * random_seq = (int *) malloc(sizeof(int)*length);
    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, length,0);
    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, length,0);
    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, length,0);
    copyIVector(initial_sequence,random_seq,length);
    copyIVector(initial_sequence,sequence_out,length);
    int improve = 1;int stay;
    //Variables for the ties
    long int it_bp=MAX_LONG;
    long int idle_time;
    int tb=0;
    int index_tie=0;
    long * array_makespan = (long *) malloc(sizeof(long)*my_flowshop.jobs);

    ftime(&actual_time);
    secs =initial_time+ actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    while(improve == 1 && secs<stopping_criterion)
    {
        improve = 0;
        //Random sequence
        index = length;
        while (index != 0)
        {
            index--;
            rnd_numb = rand()%length;
            temp = random_seq[index];
            random_seq[index] = random_seq[rnd_numb];
            random_seq[rnd_numb] = temp;
        }
        for (j = 0; j < length; j++)
        {
            copyIVector(sequence_out,sub_seq,length);
            j_aux = j;
            stay=1;
            for (k = 0; k < length && stay == 1; k++)
            {
                if (sub_seq[k] == random_seq[j])
                {
                    stay = 0;
                    j_aux = k;
                }
            }
            job_chosen = sub_seq[j_aux];
            // removes the job from
            extractIVector(sub_seq,length,j_aux);
            calculate_e_q_f(my_flowshop.machines,length-1,my_flowshop.pt,sub_seq,job_chosen,e,q, f);
            // compute best cmax
            min_local_mkspan = MAX_LONG;
            for(l=0;l<length;l++) {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++) {
                    if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < min_local_mkspan ) {
                    min_local_mkspan = best_of_position_l;
                    bp = l;
                }
            }

            // inserts the job in the best position is the solution has improved
            if(min_local_mkspan<best_makespan)
            {
                best_makespan=min_local_mkspan;
                improve=1;
                // inserts the job in the best position
                insertIVector(sub_seq,length,job_chosen,bp);
                copyIVector(sub_seq,sequence_out,length);
            }
        }
        ftime(&actual_time);
        secs =initial_time+ actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    }
    free(sub_seq);
    free(random_seq);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);
    free(array_makespan);
    return best_makespan;
}


int LS_FF(FLOWSHOP my_flowshop, int * initial_sequence,int length,int initial_makespan, float initial_time, float stopping_criterion, int * sequence_out)
{
    int i,j,k,l;
	struct timeb init_time;
	struct timeb actual_time;
	float secs;
    ftime(&init_time);
    int index,rnd_numb,temp,j_aux,job_chosen;
    int bp = -1;
    long int best_of_position_l;
    long int best_makespan = initial_makespan;
    long int min_local_mkspan=0;
    int * sub_seq = (int *) malloc(sizeof(int)*length);
    int * random_seq = (int *) malloc(sizeof(int)*length);
    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, length,0);
    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, length,0);
    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, length,0);
    copyIVector(initial_sequence,random_seq,length);
    copyIVector(initial_sequence,sequence_out,length);
    int improve = 1;int stay;
    //Variables for the ties
    long int it_bp=MAX_LONG;
    long int idle_time;
    int tb=0;
    int index_tie=0;
    long * array_makespan = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * ptb = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    long * f_prima = (long *) malloc(sizeof(long)*my_flowshop.machines);
    setval_Lvector(f_prima,my_flowshop.machines,0);

    ftime(&actual_time);
    secs =initial_time+ actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    while(improve == 1 && secs<stopping_criterion)
    {
        improve = 0;
        //Random sequence
        index = length;
        while (index != 0)
        {
            index--;
            rnd_numb = rand()%length;
            temp = random_seq[index];
            random_seq[index] = random_seq[rnd_numb];
            random_seq[rnd_numb] = temp;
        }
        for (j = 0; j < length; j++)
        {
            copyIVector(sequence_out,sub_seq,length);
            j_aux = j;
            stay=1;
            for (k = 0; k < length && stay == 1; k++)
            {
                if (sub_seq[k] == random_seq[j])
                {
                    stay = 0;
                    j_aux = k;
                }
            }
            job_chosen = sub_seq[j_aux];
            // removes the job from
            extractIVector(sub_seq,length,j_aux);
            calculate_e_q_f(my_flowshop.machines,length-1,my_flowshop.pt,sub_seq,job_chosen,e,q, f);
            // compute best cmax
            min_local_mkspan = MAX_LONG;
            for(l=0;l<length;l++) {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++) {
                    if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < min_local_mkspan ) {
                    min_local_mkspan = best_of_position_l;
                    bp = l;
                }
            }
            //number of tiesint contador = 0;
            tb=0;
            for(l=0;l<length;l++)
            {
                if (array_makespan[l] == min_local_mkspan)
                {
                    tb++;
                }
            }

            //More than 1 tie
            if(tb>1)
            {
                ptb = (int *) malloc(sizeof(int)*tb);
                index_tie=0;
                //Find the position with ties
                for(l=0;l<length;l++)
                {
                    if(array_makespan[l] == min_local_mkspan)
                    {
                        ptb[index_tie]=l;
                        index_tie++;
                    }
                }
                it_bp = MAX_LONG;
                for(l=0;l<tb;l++)
                {
                    idle_time=0;
                    if(ptb[l]<length-1)
                    {
                        f_prima[0]=f[0][ptb[l]]+(long)my_flowshop.pt[0][sub_seq[ptb[l]]];
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            f_prima[i]=lMax(f_prima[i-1],f[i][ptb[l]])+(long)my_flowshop.pt[i][sub_seq[ptb[l]]];
                            idle_time=idle_time+f[i][ptb[l]]-e[i][ptb[l]]+(long)my_flowshop.pt[i][sub_seq[ptb[l]]]+lMax(0,f_prima[i-1]-f[i][ptb[l]]);
                        }
                    }
                    else
                    {
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            idle_time=idle_time+f[i][ptb[l]]-e[i][ptb[l]-1];
                        }
                    }
                    if(idle_time<it_bp)
                    {
                        it_bp=idle_time;
                        bp=ptb[l];
                    }
                }
                // free(ptb);  // is already freed at the end of the function
            }

            // inserts the job in the best position is the solution has improved
            if(min_local_mkspan<best_makespan)
            {
                best_makespan=min_local_mkspan;
                improve=1;
                // inserts the job in the best position
                insertIVector(sub_seq,length,job_chosen,bp);
                copyIVector(sub_seq,sequence_out,length);
            }
        }
        ftime(&actual_time);
        secs =initial_time+ actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    }
    free(sub_seq);
    free(random_seq);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);
    free(array_makespan);
    free(ptb);
    free(f_prima);
    return best_makespan;
}

int IteratedGreedy_Ruiz2007_outIter(FLOWSHOP my_flowshop,int d,double T,float stopping_criterion,int* iteraciones)
{
	struct timeb init_time;
	struct timeb actual_time;
	struct timeb final_time;
	float secs;
    ftime(&init_time);
    int * sequence_ini=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_iteration=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_best_local=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_destr_constr=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int mksp_initial=NEH_Taillard_seq(my_flowshop,sequence_ini);
    int i,j,k,l;
    int random_number;double random_double,temperature_aux;
    int stay_in_loop;int number_destructed_jobs;int is_destructed;
    long int makespan_iteration = MAX_LONG;
    long int best_makespan,makespan_destr_constr,best_local_makespan;
    int N=my_flowshop.jobs;
    int M=my_flowshop.machines;
     int ** p_ij = (int **) dyn_mat(M,N,sizeof(int));
    int * best_sequence=(int *) malloc(sizeof(int)*N);
    int * destruction_jobs=(int *) malloc(sizeof(int)*N-d);
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            p_ij[i][j]=my_flowshop.pt[i][j];
        }
    }
    double suma_p_ij = 0;
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < M; i++)
        {
            suma_p_ij = suma_p_ij +(double) p_ij[i][j];
        }
    }
    double Temperature = T * suma_p_ij / (double)N / (double)M / (double)10;
    ftime(&actual_time);
    secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);

    //LOCAL SEARCH
    best_local_makespan=IterativeImprovement_Insertion(my_flowshop,sequence_ini,N,mksp_initial,secs,stopping_criterion,sequence_iteration);
    best_makespan=best_local_makespan;
    copyIVector(sequence_iteration,best_sequence,N);
    copyIVector(sequence_iteration,sequence_best_local,N);
    *iteraciones=1;
    while(secs<stopping_criterion)
    {
        *iteraciones=*iteraciones+1;
        //DESTRUCTION PHASE
        for (l = 0; l < d; l++)
        {
            destruction_jobs[l] = -1;
        }
        //Chosen random jobs
        for (l = 0; l < d; l++)
        {
            stay_in_loop = 1;
            random_number = rand()%N;
            while (stay_in_loop)
            {
                stay_in_loop = 0;
                random_number = rand()%N;
                for (i = 0; i < d && stay_in_loop == 0; i++)
                {
                    if (random_number == destruction_jobs[i])
                    {
                        stay_in_loop = 1;
                    }
                }
            }
            destruction_jobs[l] = random_number;
        }

        //CONSTRUCTION PHASE
        number_destructed_jobs = 0;
        for (j = 0; j < N; j++)
        {
            is_destructed =0;
            for (l = 0; l < d; l++)
            {
                if (sequence_iteration[j] == destruction_jobs[l])
                {
                    is_destructed = 1;
                }
            }
            if (is_destructed == 0)
            {
                sequence_destr_constr[j - number_destructed_jobs] = sequence_iteration[j];
            }
            else
            {
                number_destructed_jobs++;
            }
        }
        makespan_destr_constr=iteratedGreedy_Ruiz2007_Construction(my_flowshop,N-d,sequence_destr_constr,destruction_jobs);

        //LOCAL SEARCH
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
        best_local_makespan=IterativeImprovement_Insertion(my_flowshop,sequence_destr_constr,N,makespan_destr_constr,secs,stopping_criterion,sequence_best_local);
        random_double=(double) rand()/(double) RAND_MAX;
        if(best_local_makespan<makespan_iteration)
        {
            makespan_iteration=best_local_makespan;
            copyIVector(sequence_best_local,sequence_iteration,N);
            if(best_local_makespan<best_makespan)
            {
                best_makespan=best_local_makespan;
                copyIVector(sequence_best_local,best_sequence,N);
            }
        }
        else
        {
            temperature_aux=-(double)(best_local_makespan - makespan_iteration) / Temperature;
            if (random_double <= exp(temperature_aux))
            {
                makespan_iteration = best_local_makespan;
                copyIVector(sequence_best_local,sequence_iteration,N);
            }
        }
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    }
    free(sequence_best_local);
    free(sequence_ini);
    free(sequence_destr_constr);
    free(destruction_jobs);
    free(sequence_iteration);
    free(best_sequence);
    free_mat_int(p_ij,M);
    return best_makespan;
}

int IterativeImprovement_Insertion(FLOWSHOP my_flowshop, int* initial_sequence,int length,int initial_makespan, float initial_time, float stopping_criterion, int* sequence_out)
{
    int i,j,k,l;
	struct timeb init_time;
	struct timeb actual_time;
	float secs;
    ftime(&init_time);
    int index,rnd_numb,temp,j_aux,job_chosen,best_position;
    long int best_of_position_l;
    long int best_makespan = initial_makespan;
    long int min_local_mkspan=0;
    int * sub_seq = (int *) malloc(sizeof(int)*length);
    int * random_seq = (int *) malloc(sizeof(int)*length);
    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, length,0);
    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, length,0);
    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, length,0);
    copyIVector(initial_sequence,random_seq,length);
    copyIVector(initial_sequence,sequence_out,length);
    int improve = 1;int stay;
    ftime(&actual_time);
    secs =initial_time+ actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    while(improve == 1 && secs<stopping_criterion)
    {
        improve = 0;
        //Random sequence
        index = length;
        while (index != 0)
        {
            index--;
            rnd_numb = rand()%length;
            temp = random_seq[index];
            random_seq[index] = random_seq[rnd_numb];
            random_seq[rnd_numb] = temp;
        }
        for (j = 0; j < length; j++)
        {
            copyIVector(sequence_out,sub_seq,length);
            j_aux = j;
            stay=1;
            for (k = 0; k < length && stay == 1; k++)
            {
                if (sub_seq[k] == random_seq[j])
                {
                    stay = 0;
                    j_aux = k;
                }
            }
            job_chosen = sub_seq[j_aux];
            // removes the job from
            extractIVector(sub_seq,length,j_aux);
            calculate_e_q_f(my_flowshop.machines,length-1,my_flowshop.pt,sub_seq,job_chosen,e,q, f);
            // compute best cmax
            min_local_mkspan = MAX_LONG;
            for(l=0;l<length;l++) {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++) {
                    if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                }
                if(best_of_position_l < min_local_mkspan ) {
                    min_local_mkspan = best_of_position_l;
                    best_position = l;
                }
            }
            if(min_local_mkspan<best_makespan)
            {
                best_makespan=min_local_mkspan;
                improve=1;
                // inserts the job in the best position
                insertIVector(sub_seq,length,job_chosen,best_position);
                copyIVector(sub_seq,sequence_out,length);
            }
        }
        ftime(&actual_time);
        secs = initial_time+actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    }
    free(sub_seq);
    free(random_seq);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);
    return best_makespan;
}

int NEH_Taillard_seq(FLOWSHOP my_flowshop,int* partial_sequence) {

    register int i,j,k,l;
    long int best_of_position_l;
    long int curr_best = MAX_LONG;

    int best_position = -1;

    long * sum_processing_times = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * non_scheduled_jobs = (int *) malloc(sizeof(int)*my_flowshop.jobs);

    setval_Lvector(sum_processing_times,my_flowshop.jobs, 0);

    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);

    long * f_prima = (long *) malloc(sizeof(long)*my_flowshop.machines);
    setval_Lvector(f_prima,my_flowshop.machines,0);




    // store the sum of the completion times
    for(j=0;j<my_flowshop.jobs; j++) {
        for(i=0;i< my_flowshop.machines;i++) {
            sum_processing_times[j] += (long) my_flowshop.pt[i][j];
        }
    }

    // sort the jobs in descending order of the sum of the processing times
    brqsortLVectorD(sum_processing_times, non_scheduled_jobs,my_flowshop.jobs,'D');

    // start with the initial job (to be inserted in pos 0) and remove the job
    int tamano_actual=my_flowshop.jobs;
    partial_sequence[0] = non_scheduled_jobs[0];
    extractIVector(non_scheduled_jobs,my_flowshop.jobs,0);
    tamano_actual--;

    // for the rest of the jobs
    for(k=1;k < my_flowshop.jobs; k++) {

            // compute earliest starting time
            // first job, first machine
            e[0][0] = (long) my_flowshop.pt[0][partial_sequence[0]];
            // first job, rest of the machines
            for(i=1;i< my_flowshop.machines; i++) {
                e[i][0] = e[i-1][0] + (long) my_flowshop.pt[i][partial_sequence[0]];
            }
            // rest of the jobs
            for(j=1;j< k;j++) {
                // rest of the jobs, first machine
                e[0][j] = e[0][j-1] + (long) my_flowshop.pt[0][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=1;i< my_flowshop.machines;i++) {
                    e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute tail
            // last job, last machine
            q[my_flowshop.machines-1][k-1] = (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[k-1]];
            // last job, rest of the machines
            for(i= my_flowshop.machines-2;i>=0; i--) {
                q[i][k-1] = q[i+1][k-1] + (long) my_flowshop.pt[i][partial_sequence[k-1]];
            }
            // rest of the jobs
            for(j=k-2;j>=0;j--) {
                // rest of the jobs, last machine
                q[my_flowshop.machines-1][j] = q[my_flowshop.machines-1][j+1] + (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=my_flowshop.machines-2;i>=0;i--) {
                    q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute f (for all positions where the job can be inserted)
            // first position (l=0)
            // first position, first machine
            f[0][0] = (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
            // first position, rest of the machines
            for(i=1;i< my_flowshop.machines;i++) {
                f[i][0] = f[i-1][0] + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
            }
            // rest of the positions
            for(l=1;l<=k;l++) {
                // first machine
                f[0][l] = e[0][l-1] + (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
                // rest of the machines
                for(i=1;i<my_flowshop.machines;i++) {
                    f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
                }
            }

            // compute best cmax
            curr_best = MAX_LONG;
            for(l=0;l<=k;l++) {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++) {
                    if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                }
                if(best_of_position_l < curr_best ) {
                    curr_best = best_of_position_l;
                    best_position = l;
                }
            }

            // inserts the job in the best position
            insertIVector(partial_sequence,k+1,non_scheduled_jobs[0],best_position);

            // removes the job from the non_scheduled
            extractIVector(non_scheduled_jobs,tamano_actual,0);
            tamano_actual--;

    }


    free(non_scheduled_jobs);
    free(sum_processing_times);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);

    return curr_best;
}

int iteratedGreedy_Ruiz2007_Construction(FLOWSHOP my_flowshop,int initial_length,int* partial_sequence,int* destructed_jobs)
{

    register int i,j,k,l;
    long int best_of_position_l;
    long int curr_best = MAX_LONG;

    int best_position = -1;

    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);


    // for the rest of the jobs
    for(k=initial_length;k < my_flowshop.jobs; k++) {

            // compute earliest starting time
            // first job, first machine
            e[0][0] = (long) my_flowshop.pt[0][partial_sequence[0]];
            // first job, rest of the machines
            for(i=1;i< my_flowshop.machines; i++) {
                e[i][0] = e[i-1][0] + (long) my_flowshop.pt[i][partial_sequence[0]];
            }
            // rest of the jobs
            for(j=1;j< k;j++) {
                // rest of the jobs, first machine
                e[0][j] = e[0][j-1] + (long) my_flowshop.pt[0][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=1;i< my_flowshop.machines;i++) {
                    e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute tail
            // last job, last machine
            q[my_flowshop.machines-1][k-1] = (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[k-1]];
            // last job, rest of the machines
            for(i= my_flowshop.machines-2;i>=0; i--) {
                q[i][k-1] = q[i+1][k-1] + (long) my_flowshop.pt[i][partial_sequence[k-1]];
            }
            // rest of the jobs
            for(j=k-2;j>=0;j--) {
                // rest of the jobs, last machine
                q[my_flowshop.machines-1][j] = q[my_flowshop.machines-1][j+1] + (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[j]];
                // rest of the jobs, rest of the machines
                for(i=my_flowshop.machines-2;i>=0;i--) {
                    q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
                }
            }

            // compute f (for all positions where the job can be inserted)
            // first position (l=0)
            // first position, first machine
            f[0][0] = (long) my_flowshop.pt[0][destructed_jobs[k-initial_length]];
            // first position, rest of the machines
            for(i=1;i< my_flowshop.machines;i++) {
                f[i][0] = f[i-1][0] + (long) my_flowshop.pt[i][destructed_jobs[k-initial_length]];
            }
            // rest of the positions
            for(l=1;l<=k;l++) {
                // first machine
                f[0][l] = e[0][l-1] + (long) my_flowshop.pt[0][destructed_jobs[k-initial_length]];
                // rest of the machines
                for(i=1;i<my_flowshop.machines;i++) {
                    f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) my_flowshop.pt[i][destructed_jobs[k-initial_length]];
                }
            }

            // compute best cmax
            curr_best = MAX_LONG;
            for(l=0;l<=k;l++) {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++) {
                    if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                }
                if(best_of_position_l < curr_best ) {
                    curr_best = best_of_position_l;
                    best_position = l;
                }
            }

            // inserts the job in the best position
            insertIVector(partial_sequence,k+1,destructed_jobs[k-initial_length],best_position);


    }


    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);

    return curr_best;
}

int IG_BoB(FLOWSHOP my_flowshop,int d,double T,float stopping_criterion,int* iteraciones)
{
	struct timeb init_time;
	struct timeb actual_time;
	float secs;
    ftime(&init_time);
    int * sequence_ini=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_iteration=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_best_local=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sequence_destr_constr=(int *) malloc(sizeof(int)*my_flowshop.jobs);

    int * sequence_destr_constr_LS=(int *) malloc(sizeof(int)*my_flowshop.jobs);
    int mksp_initial=BR_R_FF(my_flowshop,sequence_ini);
    int i,j,l;
    int random_number;double random_double,temperature_aux;
    int stay_in_loop;int number_destructed_jobs;int is_destructed;
    long int best_makespan, makespan_destr_constr,best_local_makespan;
    long int makespan_iteration = MAX_LONG;
    int N=my_flowshop.jobs;
    int M=my_flowshop.machines;
    long int ** p_ij = (long int **) dyn_mat(M,N,sizeof(long int));
    int * best_sequence=(int *) malloc(sizeof(int)*N);
    int * destruction_jobs=(int *) malloc(sizeof(int)*N);
    for(i=0;i<M;i++)
    {
        for(j=0;j<N;j++)
        {
            p_ij[i][j]=my_flowshop.pt[i][j];
        }
    }
    double suma_p_ij = 0;
    for (j = 0; j < N; j++)
    {
        for (i = 0; i < M; i++)
        {
            suma_p_ij = suma_p_ij +(double) p_ij[i][j];
        }
    }
    double Temperature = T * suma_p_ij / (double)N / (double)M / (double)10;
    ftime(&actual_time);
    secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);


    //LOCAL SEARCH
    best_local_makespan=Pc_FF(my_flowshop,sequence_ini,N,mksp_initial,secs,stopping_criterion,sequence_iteration);

    best_makespan=best_local_makespan;
    copyIVector(sequence_iteration,best_sequence,N);
    copyIVector(sequence_iteration,sequence_best_local,N);
    *iteraciones=0;
    while(secs<stopping_criterion)
    {
        *iteraciones=*iteraciones+1;
        //DESTRUCTION PHASE
        for (l = 0; l < d; l++)
        {
            destruction_jobs[l] = -1;
        }
        //Chosen random jobs
        for (l = 0; l < d; l++)
        {
            stay_in_loop = 1;
            random_number = rand()%N;
            while (stay_in_loop)
            {
                stay_in_loop = 0;
                random_number = rand()%N;
                for (i = 0; i < d && stay_in_loop == 0; i++)
                {
                    if (random_number == destruction_jobs[i])
                    {
                        stay_in_loop = 1;
                    }
                }
            }
            destruction_jobs[l] = random_number;
        }

        //CONSTRUCTION PHASE
        number_destructed_jobs = 0;
        for (j = 0; j < N; j++)
        {
            is_destructed =0;
            for (l = 0; l < d; l++)
            {
                if (sequence_iteration[j] == destruction_jobs[l])
                {
                    is_destructed = 1;
                }
            }
            if (is_destructed == 0)
            {
                sequence_destr_constr[j - number_destructed_jobs] = sequence_iteration[j];
            }
            else
            {
                number_destructed_jobs++;
            }
        }
        makespan_destr_constr=Makespan_withoutTaillard_PFSP_part(my_flowshop,N-d,sequence_destr_constr);
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
        //LOCAL SEARCH
        LS_FF(my_flowshop,sequence_destr_constr,N-d,makespan_destr_constr,secs,stopping_criterion,sequence_destr_constr_LS);
        //CONSTRUCTION
        makespan_destr_constr=Construction_FF(my_flowshop,N-d,sequence_destr_constr_LS,destruction_jobs);
        //LOCAL SEARCH
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
        best_local_makespan=Pc_FF(my_flowshop,sequence_destr_constr_LS,N,makespan_destr_constr,secs,stopping_criterion,sequence_best_local);
        random_double=(double) rand()/(double) RAND_MAX;
        // printf("bestlocal:%d\t best:%d\n", best_local_makespan, makespan_iteration);
        if(best_local_makespan<makespan_iteration)
        {
            makespan_iteration=best_local_makespan;
            copyIVector(sequence_best_local,sequence_iteration,N);
            if(best_local_makespan<best_makespan)
            {
                best_makespan=best_local_makespan;
                copyIVector(sequence_best_local,best_sequence,N);
            }
        }
        else
        {
            temperature_aux=-(double)(best_local_makespan - makespan_iteration) / Temperature;
            if (random_double <= exp(temperature_aux))
            {
                makespan_iteration = best_local_makespan;
                copyIVector(sequence_best_local,sequence_iteration,N);
            }
        }
        ftime(&actual_time);
        secs = actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    }
    free(sequence_best_local);
    free(sequence_ini);
    free(sequence_destr_constr);
    free(destruction_jobs);
    free(sequence_iteration);
    free(best_sequence);
    free(sequence_destr_constr_LS);
    free_mat_int(p_ij,M);
    return best_makespan;
}

int Pc_FF(FLOWSHOP my_flowshop, int * initial_sequence,int length,int initial_makespan, float initial_time, float stopping_criterion, int * sequence_out)
{
    int i,j,k,l;
	struct timeb init_time;
	struct timeb actual_time;
	float secs;
    ftime(&init_time);
    long int index,rnd_numb,temp,j_aux,job_chosen1,job_chosen2,n_elems,number_ties,random_number;
    int bp = -1;
    long int best_of_position_l;
    long int best_makespan = initial_makespan;
    long int min_local_mkspan=0;
    int * sub_seq = (int *) malloc(sizeof(int)*length);
    long int lengthC=my_flowshop.machines*length*2;
    int * random_seq = (int *) malloc(sizeof(int)*lengthC);
    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);
    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);
    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);
    int ** trab_selecc= (int **) dyn_mat(lengthC,2, sizeof(int));
    setval_Lmatrix(trab_selecc,lengthC, 2,0);
    copyIVector(initial_sequence,sequence_out,length);
    int improve = 1;int stay=1;
    //Variables for the ties
    long int it_bp=MAX_LONG;
    long int idle_time;
    int index_tie=0;
    long int * array_makespan = (long int *) malloc(sizeof(long int)*my_flowshop.jobs);

    int * ptb = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    long * f_prima = (long *) malloc(sizeof(long)*my_flowshop.machines);
    setval_Lvector(f_prima,my_flowshop.machines,0);
    ftime(&actual_time);
    secs =initial_time+ actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    while(improve == 1 && secs<stopping_criterion)
    {
        improve = 0;
        //Critical path
        n_elems=calculate_e_q(my_flowshop.machines,length,my_flowshop.pt,sequence_out,e,q,trab_selecc);
        //Random sequence
        for (j = 0; j < n_elems; j++)
        {
            random_seq[j]=j;
        }
        index = n_elems;
        while (index != 0)
        {
            index--;
            rnd_numb = rand()%n_elems;
            temp = random_seq[index];
            random_seq[index] = random_seq[rnd_numb];
            random_seq[rnd_numb] = temp;
        }
        for (j = 0; j < n_elems; j++)
        {
            copyIVector(sequence_out,sub_seq,length);
            // removes first job from
            j_aux = j;
            stay=1;
            for (k = 0; k < length && stay == 1; k++)
            {
                if (sub_seq[k] == trab_selecc[random_seq[j]][0])
                {
                    stay = 0;
                    j_aux = k;
                }
            }
            job_chosen1 = sub_seq[j_aux];
            extractIVector(sub_seq,length,j_aux);
            // removes second job from
            j_aux = 0;
            stay=1;
            for (k = 0; k < length-1 && stay == 1; k++)
            {
                if (sub_seq[k] == trab_selecc[random_seq[j]][1])
                {
                    stay = 0;
                    j_aux = k;
                }
            }
            job_chosen2 = sub_seq[j_aux];
            extractIVector(sub_seq,length-1,j_aux);
            //Insert first job
            calculate_e_q_f(my_flowshop.machines,length-2,my_flowshop.pt,sub_seq,job_chosen1,e,q, f);
            min_local_mkspan = MAX_LONG;
            for(l=0;l<length-1;l++) {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++) {
                    if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < min_local_mkspan ) {
                    min_local_mkspan = best_of_position_l;
                    bp = l;
                }
            }
            number_ties=0;
            for(l=0;l<length-1;l++)
            {
                if (array_makespan[l] == min_local_mkspan)
                {
                    number_ties++;
                }
            }
            if(number_ties>1)
            {
                index_tie=0;
                //Find the position with ties
                for(l=0;l<length-1;l++)
                {
                    if(array_makespan[l] == min_local_mkspan)
                    {
                        ptb[index_tie]=l;
                        index_tie++;
                    }
                }
                it_bp = MAX_LONG;
                for(l=0;l<number_ties;l++)
                {
                    idle_time=0;
                    if(ptb[l]<length-2)
                    {
                        f_prima[0]=f[0][ptb[l]]+(long)my_flowshop.pt[0][sub_seq[ptb[l]]];
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            f_prima[i]=lMax(f_prima[i-1],f[i][ptb[l]])+(long)my_flowshop.pt[i][sub_seq[ptb[l]]];
                            idle_time=idle_time+f[i][ptb[l]]-e[i][ptb[l]]+(long)my_flowshop.pt[i][sub_seq[ptb[l]]]+lMax(0,f_prima[i-1]-f[i][ptb[l]]);
                        }
                    }
                    else
                    {
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            idle_time=idle_time+f[i][ptb[l]]-e[i][ptb[l]-1];
                        }
                    }
                    if(idle_time<it_bp)
                    {
                        it_bp=idle_time;
                        bp=ptb[l];
                    }
                }
            }
            insertIVector(sub_seq,length-1,job_chosen1,bp);
            //Insert second job
            calculate_e_q_f(my_flowshop.machines,length-1,my_flowshop.pt,sub_seq,job_chosen2,e,q, f);
            min_local_mkspan = MAX_LONG;
            for(l=0;l<length;l++) {
                best_of_position_l = 0;
                for(i=0;i< my_flowshop.machines;i++) {
                    if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                }
                array_makespan[l]=best_of_position_l;
                if(best_of_position_l < min_local_mkspan ) {
                    min_local_mkspan = best_of_position_l;
                    bp = l;
                }
            }
            number_ties=0;
            for(l=0;l<length;l++)
            {
                if (array_makespan[l] == min_local_mkspan)
                {
                    number_ties++;
                }
            }
            if(number_ties>1)
            {
                index_tie=0;
                //Find the position with ties
                for(l=0;l<length;l++)
                {
                    if(array_makespan[l] == min_local_mkspan)
                    {
                        ptb[index_tie]=l;
                        index_tie++;
                    }
                }
                it_bp = MAX_LONG;
                for(l=0;l<number_ties;l++)
                {
                    idle_time=0;
                    if(ptb[l]<length-1)
                    {
                        f_prima[0]=f[0][ptb[l]]+(long)my_flowshop.pt[0][sub_seq[ptb[l]]];
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            f_prima[i]=lMax(f_prima[i-1],f[i][ptb[l]])+(long)my_flowshop.pt[i][sub_seq[ptb[l]]];
                            idle_time=idle_time+f[i][ptb[l]]-e[i][ptb[l]]+(long)my_flowshop.pt[i][sub_seq[ptb[l]]]+lMax(0,f_prima[i-1]-f[i][ptb[l]]);
                        }
                    }
                    else
                    {
                        for(i=1;i<my_flowshop.machines;i++)
                        {
                            idle_time=idle_time+f[i][ptb[l]]-e[i][ptb[l]-1];
                        }
                    }
                    if(idle_time<it_bp)
                    {
                        it_bp=idle_time;
                        bp=ptb[l];
                    }
                }
            }
            // inserts the job in the best position is the solution has improved
            if(min_local_mkspan<best_makespan)
            {
                best_makespan=min_local_mkspan;
                improve=1;
                // inserts the job in the best position
                insertIVector(sub_seq,length,job_chosen2,bp);
                copyIVector(sub_seq,sequence_out,length);
            }
        }
        ftime(&actual_time);
        secs =initial_time+ actual_time.time - init_time.time + (((float) (actual_time.millitm - init_time.millitm))/1000);
    }
    free(f_prima);
    free(ptb);
    free(sub_seq);
    free(random_seq);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);
    free_mat_long(trab_selecc, lengthC);
    free(array_makespan);
    return best_makespan;
}

int Makespan_withoutTaillard_PFSP_part(FLOWSHOP my_flowshop,int n_aux, int* seq)
{
    int Makespan = 0;
    long int * load_m = (long int *) malloc(sizeof(long int)*my_flowshop.machines);
    int i,j;
    for ( i = 0; i < my_flowshop.machines; i++)
    {
        load_m[i] = 0;
    }
    for ( j = 0; j < n_aux; j++)
    {
        for ( i = 0; i < my_flowshop.machines; i++)
        {
            if (i == 0)
            {
                load_m[i] = load_m[i] + my_flowshop.pt[i][ seq[j]];
            }
            else
            {
                if (load_m[i - 1] >= load_m[i])
                {
                    load_m[i] = load_m[i - 1] + my_flowshop.pt[i][ seq[j]];
                }
                else
                {
                    load_m[i] = load_m[i] + my_flowshop.pt[i][ seq[j]];
                }
            }
        }
    }
    Makespan = load_m[my_flowshop.machines - 1];
    free(load_m);
    return Makespan;
}

int calculate_e_q(int m,int k, int ** p_ij,int * partial_sequence,long int ** e,long int ** q,int ** trab_selecc)
{
    int i,j,l;
    // compute earliest starting time
    // first job, first machine
    e[0][0] = (long) p_ij[0][partial_sequence[0]];
    // first job, rest of the machines
    for(i=1;i< m; i++) {
        e[i][0] = e[i-1][0] + (long) p_ij[i][partial_sequence[0]];
    }
    // rest of the jobs
    for(j=1;j< k;j++) {
        // rest of the jobs, first machine
        e[0][j] = e[0][j-1] + (long) p_ij[0][partial_sequence[j]];
        // rest of the jobs, rest of the machines
        for(i=1;i< m;i++) {
            e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) p_ij[i][partial_sequence[j]];
        }
    }
    // compute tail
    // last job, last machine
    q[m-1][k-1] = (long) p_ij[m-1][partial_sequence[k-1]];
    // last job, rest of the machines
    for(i= m-2;i>=0; i--) {
        q[i][k-1] = q[i+1][k-1] + (long) p_ij[i][partial_sequence[k-1]];
    }
    // rest of the jobs
    for(j=k-2;j>=0;j--) {
        // rest of the jobs, last machine
        q[m-1][j] = q[m-1][j+1] + (long) p_ij[m-1][partial_sequence[j]];
        // rest of the jobs, rest of the machines
        for(i=m-2;i>=0;i--) {
            q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) p_ij[i][partial_sequence[j]];
        }
    }
    int Cmax=e[m-1][k-1];
    int elements=0;
    for(j=0;j<k; j++)
    {
        for(i=0;i< m-1; i++)
        {
            if(e[i][j]+q[i+1][j]==Cmax)
            {
                if(j>0)
                {
                    trab_selecc[elements][0]=partial_sequence[j-1];
                    trab_selecc[elements][1]=partial_sequence[j];
                    elements++;
                }
                if(j<k-1)
                {
                    trab_selecc[elements][0]=partial_sequence[j];
                    trab_selecc[elements][1]=partial_sequence[j+1];
                    elements++;
                }
                break;
            }
        }
    }
    return elements;
}


int BR_R_FF(FLOWSHOP my_flowshop,int* partial_sequence) {

    register int i,j,k,l;
    long int best_of_position_l;
    long int curr_best = MAX_LONG;
    long int it_bp=MAX_LONG;
    long int idle_time;

    int bp = -1,number_ties,random_number,temp,ind;
    int iii,length,job_chosen1,job_chosen2,min_local_mkspan,j_aux,stay;
    int index=0;

    long * sum_processing_times = (long *) malloc(sizeof(long)*my_flowshop.jobs);
    int * non_scheduled_jobs = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    int * sub_seq = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    long * array_makespan = (long *) malloc(sizeof(long)*my_flowshop.jobs);

    setval_Lvector(sum_processing_times,my_flowshop.jobs, 0);

    long int ** e = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(e,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** q = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(q,my_flowshop.machines, my_flowshop.jobs,0);

    long int ** f = (long int **) dyn_mat(my_flowshop.machines, my_flowshop.jobs, sizeof(long int));
    setval_Lmatrix(f,my_flowshop.machines, my_flowshop.jobs,0);

    // store the sum of the completion times
    for(j=0;j<my_flowshop.jobs; j++) {
        for(i=0;i< my_flowshop.machines;i++) {
            sum_processing_times[j] += (long) my_flowshop.pt[i][j];
        }
    }

    // sort the jobs in descending order of the sum of the processing times
    for(j=0;j< my_flowshop.jobs;j++) non_scheduled_jobs[j] = j;
    brqsortLVectorD(sum_processing_times, non_scheduled_jobs, 0, my_flowshop.jobs-1);


    int * ptb = (int *) malloc(sizeof(int)*my_flowshop.jobs);
    long * f_prima = (long *) malloc(sizeof(long)*my_flowshop.machines);
    setval_Lvector(f_prima,my_flowshop.machines,0);


    // start with the initial job (to be inserted in pos 0) and remove the job
    int tamano_actual=my_flowshop.jobs;
    partial_sequence[0] = non_scheduled_jobs[0];
    extractIVector(non_scheduled_jobs,my_flowshop.jobs,0);
    tamano_actual--;

    // for the rest of the jobs
    for(k=1;k < my_flowshop.jobs; k++) {

        // compute earliest starting time
        // first job, first machine
        e[0][0] = (long) my_flowshop.pt[0][partial_sequence[0]];
        // first job, rest of the machines
        for(i=1;i< my_flowshop.machines; i++) {
            e[i][0] = e[i-1][0] + (long) my_flowshop.pt[i][partial_sequence[0]];
        }
        // rest of the jobs
        for(j=1;j< k;j++) {
            // rest of the jobs, first machine
            e[0][j] = e[0][j-1] + (long) my_flowshop.pt[0][partial_sequence[j]];
            // rest of the jobs, rest of the machines
            for(i=1;i< my_flowshop.machines;i++) {
                e[i][j] = lMax(e[i-1][j], e[i][j-1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
            }
        }

        // compute tail
        // last job, last machine
        q[my_flowshop.machines-1][k-1] = (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[k-1]];
        // last job, rest of the machines
        for(i= my_flowshop.machines-2;i>=0; i--) {
            q[i][k-1] = q[i+1][k-1] + (long) my_flowshop.pt[i][partial_sequence[k-1]];
        }
        // rest of the jobs
        for(j=k-2;j>=0;j--) {
            // rest of the jobs, last machine
            q[my_flowshop.machines-1][j] = q[my_flowshop.machines-1][j+1] + (long) my_flowshop.pt[my_flowshop.machines-1][partial_sequence[j]];
            // rest of the jobs, rest of the machines
            for(i=my_flowshop.machines-2;i>=0;i--) {
                q[i][j] = lMax(q[i+1][j], q[i][j+1]) + (long) my_flowshop.pt[i][partial_sequence[j]];
            }
        }

        // compute f (for all positions where the job can be inserted)
        // first position (l=0)
        // first position, first machine
        f[0][0] = (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
        // first position, rest of the machines
        for(i=1;i< my_flowshop.machines;i++) {
            f[i][0] = f[i-1][0] + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
        }
        // rest of the positions
        for(l=1;l<=k;l++) {
            // first machine
            f[0][l] = e[0][l-1] + (long) my_flowshop.pt[0][non_scheduled_jobs[0]];
            // rest of the machines
            for(i=1;i<my_flowshop.machines;i++) {
                f[i][l] = lMax(e[i][l-1], f[i-1][l]) + (long) my_flowshop.pt[i][non_scheduled_jobs[0]];
            }
        }

        // compute best cmax
        curr_best = MAX_LONG;
        for(l=0;l<=k;l++)
        {
            best_of_position_l = 0;
            for(i=0;i< my_flowshop.machines;i++)
            {
                if( f[i][l] + q[i][l] > best_of_position_l)
                best_of_position_l = f[i][l] + q[i][l] ;
            }
            array_makespan[l]=best_of_position_l;
            if(best_of_position_l < curr_best )
            {
                curr_best = best_of_position_l;
                bp = l;
            }
        }
        //number of ties
        number_ties=0;
        for(l=0;l<=k;l++)
        {
            if (array_makespan[l] == curr_best)
            {
                number_ties++;
            }
        }

        //More than 1 tie
        if(number_ties>1)
        {
            index=0;
            //Find the position with ties
            for(l=0;l<=k;l++)
            {
                if(array_makespan[l] == curr_best)
                {
                    ptb[index]=l;
                    index++;
                }
            }
            it_bp = MAX_LONG;
            for(j=0;j<number_ties;j++)
            {
                idle_time=0;
                if(ptb[j]<k)
                {
                    f_prima[0]=f[0][ptb[j]]+(long)my_flowshop.pt[0][partial_sequence[ptb[j]]];
                    for(i=1;i<my_flowshop.machines;i++)
                    {
                        f_prima[i]=lMax(f_prima[i-1],f[i][ptb[j]])+(long)my_flowshop.pt[i][partial_sequence[ptb[j]]];
                        idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]]+(long)my_flowshop.pt[i][partial_sequence[ptb[j]]]+lMax(0,f_prima[i-1]-f[i][ptb[j]]);
                    }
                }
                else
                {
                    for(i=1;i<my_flowshop.machines;i++)
                    {
                        idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]-1];
                    }
                }
                if(idle_time<it_bp)
                {
                    it_bp=idle_time;
                    bp=ptb[j];
                }
            }
        }

        // inserts the job in the best position
        insertIVector(partial_sequence,k+1,non_scheduled_jobs[0],bp);
        length=k+1;
        // removes the job from the non_scheduled
        extractIVector(non_scheduled_jobs,tamano_actual,0);
        if(k>3)
        {
            if(((k+1)/2-1)*2+1<k+1)
            {
                for(iii=0;iii<(k+1)/2;iii++)
                {
                    //int numAleat=rand()%2;
                    copyIVector(partial_sequence,sub_seq,length);
                    // removes first job from
                    j_aux = iii*2;
                    stay=1;
                    for (ind= 0; ind< length && stay == 1; ind++)
                    {
                        if (sub_seq[ind] == partial_sequence[iii*2])//if (sub_seq[k] == trab_selecc[random_seq[j]][numAleat])
                        {
                            stay = 0;
                            j_aux = ind;
                        }
                    }
                    job_chosen1 = sub_seq[j_aux];
                    extractIVector(sub_seq,length,j_aux);
                    // removes second job from
                    stay=1;
                    for (ind = 0; ind < length-1 && stay == 1; ind++)
                    {
                        if (sub_seq[ind] == partial_sequence[iii*2+1])//if (sub_seq[k] == trab_selecc[random_seq[j]][1-numAleat])
                        {
                            stay = 0;
                            j_aux = ind;
                        }
                    }
                    job_chosen2 = sub_seq[j_aux];
                    extractIVector(sub_seq,length-1,j_aux);
                    //Insert first job
                    calculate_e_q_f(my_flowshop.machines,length-2,my_flowshop.pt,sub_seq,job_chosen1,e,q, f);
                    min_local_mkspan = MAX_LONG;
                    for(l=0;l<length-1;l++) {
                        best_of_position_l = 0;
                        for(i=0;i< my_flowshop.machines;i++) {
                            if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                        }
                        array_makespan[l]=best_of_position_l;
                        if(best_of_position_l < min_local_mkspan ) {
                            min_local_mkspan = best_of_position_l;
                            bp = l;
                        }
                    }
                    number_ties=0;
                    for(l=0;l<length-1;l++)
                    {
                        if (array_makespan[l] == min_local_mkspan)
                        {
                            number_ties++;
                        }
                    }
                    if(number_ties>1)
                    {
                        index=0;
                        //Find the position with ties
                        for(l=0;l<length-1;l++)
                        {
                            if(array_makespan[l] == min_local_mkspan)
                            {
                                ptb[index]=l;
                                index++;
                            }
                        }
                        it_bp = MAX_LONG;
                        for(j=0;j<number_ties;j++)
                        {
                            idle_time=0;
                            if(ptb[j]<length-2)
                            {
                                f_prima[0]=f[0][ptb[j]]+(long)my_flowshop.pt[0][sub_seq[ptb[j]]];
                                for(i=1;i<my_flowshop.machines;i++)
                                {
                                    f_prima[i]=lMax(f_prima[i-1],f[i][ptb[j]])+(long)my_flowshop.pt[i][sub_seq[ptb[j]]];
                                    idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]]+(long)my_flowshop.pt[i][sub_seq[ptb[j]]]+lMax(0,f_prima[i-1]-f[i][ptb[j]]);
                                }
                            }
                            else
                            {
                                for(i=1;i<my_flowshop.machines;i++)
                                {
                                    idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]-1];
                                }
                            }
                            if(idle_time<it_bp)
                            {
                                it_bp=idle_time;
                                bp=ptb[j];
                            }
                        }
                    }
                    insertIVector(sub_seq,length-1,job_chosen1,bp);
                    //Insert second job
                    calculate_e_q_f(my_flowshop.machines,length-1,my_flowshop.pt,sub_seq,job_chosen2,e,q, f);
                    min_local_mkspan = MAX_LONG;
                    for(l=0;l<length;l++) {
                        best_of_position_l = 0;
                        for(i=0;i< my_flowshop.machines;i++) {
                            if( f[i][l] + q[i][l] > best_of_position_l) best_of_position_l = f[i][l] + q[i][l] ;
                        }
                        array_makespan[l]=best_of_position_l;
                        if(best_of_position_l < min_local_mkspan ) {
                            min_local_mkspan = best_of_position_l;
                            bp = l;
                        }
                    }
                    number_ties=0;
                    for(l=0;l<length;l++)
                    {
                        if (array_makespan[l] == min_local_mkspan)
                        {
                            number_ties++;
                        }
                    }
                    if(number_ties>1)
                    {
                        index=0;
                        //Find the position with ties
                        for(l=0;l<length;l++)
                        {
                            if(array_makespan[l] == min_local_mkspan)
                            {
                                ptb[index]=l;
                                index++;
                            }
                        }
                        it_bp = MAX_LONG;
                        for(j=0;j<number_ties;j++)
                        {
                            idle_time=0;
                            if(ptb[j]<length-1)
                            {
                                f_prima[0]=f[0][ptb[j]]+(long)my_flowshop.pt[0][sub_seq[ptb[j]]];
                                for(i=1;i<my_flowshop.machines;i++)
                                {
                                    f_prima[i]=lMax(f_prima[i-1],f[i][ptb[j]])+(long)my_flowshop.pt[i][sub_seq[ptb[j]]];
                                    idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]]+(long)my_flowshop.pt[i][sub_seq[ptb[j]]]+lMax(0,f_prima[i-1]-f[i][ptb[j]]);
                                }
                            }
                            else
                            {
                                for(i=1;i<my_flowshop.machines;i++)
                                {
                                    idle_time=idle_time+f[i][ptb[j]]-e[i][ptb[j]-1];
                                }
                            }
                            if(idle_time<it_bp)
                            {
                                it_bp=idle_time;
                                bp=ptb[j];
                            }
                        }
                    }
                    // inserts the job in the best position is the solution has improved
                    if(min_local_mkspan<curr_best)
                    {
                        curr_best=min_local_mkspan;
                        // inserts the job in the best position
                        insertIVector(sub_seq,length,job_chosen2,bp);
                        copyIVector(sub_seq,partial_sequence,length);
                    }
                }
                tamano_actual--;
            }
        }
    }
    free(non_scheduled_jobs);
    free(sum_processing_times);
    free(array_makespan);
    free(sub_seq);
    free_mat_long(e, my_flowshop.machines);
    free_mat_long(q, my_flowshop.machines);
    free_mat_long(f, my_flowshop.machines);

    return curr_best;
}



int main(int argc, char *argv[]) {
    if ( argc < 4 ) {
        printf("USAGE: %s INST T OUTPUT.json\n", argv[0]);
        printf("\tT: timelimit parameter (computingTime = T * nbjobs * nbmachines / 2 milliseconds)\n");
        exit(1);
    }
    srand(0);  // set seed to 0
    double t = atof(argv[2]);
	struct timeb init_time;
	struct timeb final_time;
	float secs;
    int i,j;
    FLOWSHOP my_flowshop;
    int d;double T;
    const int runs=10;
    my_flowshop.pt = loadPTimes_nrows(argv[1], &my_flowshop.jobs, &my_flowshop.machines);
    int length=3;
    int objective_array[runs];
    double AMK = 0.;
    double ACPU = 0.;
    int it1=0;
    double iterations1=0;
    double stopping_criterion = (double)my_flowshop.jobs * (double)my_flowshop.machines/ (double)2 * (double)t / (double)1000;
    printf("### IGBOB: running experiments for: (%s,%d)\n", argv[1], (int)t);
    printf("computation time %d*n*m/2 = %f\n", (int) t, stopping_criterion);
    int obj_array[runs];
    for(i=0;i<runs;i++)
    {
        clock_t instant_time;
        j=0;
        int MK_IG=0;

        //Iterated greedy Benavides
        it1=0;
        ftime(&init_time);
        d=2;
        T=0.7;
        MK_IG=IG_BoB(my_flowshop,d,T,stopping_criterion,&it1);
        obj_array[i] = (int)MK_IG;
        iterations1+=it1;
        ftime(&final_time);
        secs = final_time.time - init_time.time + (((float) (final_time.millitm - init_time.millitm))/1000);
        AMK += (double)MK_IG;
        ACPU += secs;
        printf("run %d objective:\t%d\n",i,MK_IG);
        j++;
    }
    printf("END IG_BoB\n");
    printf("avg makespan:\t%f\n",(AMK/(double)runs));
    printf("avg iterations:\t%f\n",iterations1/(double)runs);
    // write experiment recap
    FILE* f = fopen(argv[3], "w");
    fprintf(f, "{\n");
    fprintf(f, "\t\"inst_name\":\"%s\",\n",argv[1]);
    fprintf(f, "\t\"time_searched\":%f,\n",stopping_criterion);
    fprintf(f, "\t\"average_primal\":\"%f\",\n",(AMK/(double)runs));
    fprintf(f, "\t\"average_nb_iters\":\"%f\",\n",iterations1/(double)runs);
    fprintf(f, "\t\"primal_list\": [");
        for ( i = 0 ; i < runs ; i++ ) {
            if ( i < runs-1 ) {
                fprintf(f, "%d, ", obj_array[i]);
            } else {
                fprintf(f, "%d", obj_array[i]);
            }
        }
    fprintf(f, "]\n");
    fprintf(f, "}");
    fclose(f);
    return 0;
}
