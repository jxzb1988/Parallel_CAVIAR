#include <vector>
#include <algorithm>
#include <set>
#include <iostream>
#include <armadillo>
#include <omp.h>
#include "Util.h"
#include "PostCal.h"
//#include <array>
using namespace arma;


void printGSLPrint(mat &A, int row, int col) {
	for(int i = 0; i < row; i++) {
		for(int j = 0; j < col; j++)
			printf("%g ", A(i, j));
		printf("\n");
	}	
}

string PostCal::convertConfig2String(int * config, int size) {
	string result = "0";
	for(int i = 0; i < size; i++)
		if(config[i]==1)
			result+= "_" + convertInt(i);
	return result;
}

//double PostCal::likelihood(int * configure, double * stat, double NCP) {
double PostCal::likelihood(int * configure_x, double NCP) {



    //    cout<<"Calculate likelihood of a configure"<<endl;
    //    cout<<"Value of NCP is: "<<NCP<<endl;
	int causalCount = 0;
	int index_C = 0;
        double matDet = 0;
	double res    = 0;

	for(int i = 0; i < snpCount; i++)
         { 
           causalCount += configure_x[i];
       //    cout<<" "<<configure_x[i];
         }
      //  cout<<endl;
      //  cout<<"Value of causalCount is: "<<causalCount<<endl;
       // return (0.1);
	if(causalCount == 0){
		mat tmpResultMatrix1N = statMatrixtTran * invSigmaMatrix;
		mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
		res = tmpResultMatrix11(0,0);	
		baseValue = res;
		matDet = sigmaDet;
		res = res - baseValue;
		return( exp(-res/2)/sqrt(abs(matDet)) );
	}
      //  return (0.1);
	mat U(snpCount, causalCount, fill::zeros);
	mat V(causalCount, snpCount, fill::zeros);
	mat VU(causalCount, causalCount, fill::zeros);
      //	return (0.1);

	for(int i = 0; i < snpCount; i++) {
                if (configure_x[i] == 0)	continue;
                else {
//                        cout<<"Value of index_C is: "<<index_C<<endl;
                        if(index_C+1>causalCount)
                         {
         //                  cout<<"Value of index_C is: "<<index_C<<" and value of causalCount is: "<<causalCount<<endl;
        //                   cout<<"Error, return"<<endl;
        //                   return (0.1);
                         }
                        for(int j = 0; j < snpCount; j++){ 
                              
                                U(j, index_C) = sigmaMatrix(j,i);
                        }
                      //          U(j, 0) = sigmaMatrix(j,i);
//                          {
//                             cout<<"Output value: "<<index_C<<" "<<j<<endl;
//                             U(j, index_C) = 2;
//                          }
			V(index_C, i) = NCP;
           //             cout<<"Value of index_C is: "<<index_C<<endl;
                        if(configure_x[i] == 1)
                         {
                           index_C++;
                         }
                }
        }
      //  return (0.1);
	VU = V * U;
     //   return (0.1);
	mat I_AA   = mat(snpCount, snpCount, fill::eye);
	mat tmp_CC = mat(causalCount, causalCount, fill::eye)+ VU;
	matDet = det(tmp_CC) * sigmaDet;
	mat tmp_AA = invSigmaMatrix - (invSigmaMatrix * U) * pinv(tmp_CC) * V ;
     //   return (0.1); 
	//tmp_AA     = invSigmaMatrix * tmp_AA;
	mat tmpResultMatrix1N = statMatrixtTran * tmp_AA;
        mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
        res = tmpResultMatrix11(0,0);  

	res = res - baseValue;
	if(matDet==0) {
//		cout << "Error the matrix is singular and we fail to fix it." << endl;
		exit(0);
	}
	/*
		We compute the log of -res/2-log(det) to see if it is too big or not. 
		In the case it is too big we just make it a MAX value.
	*/
     //   delete(U);
     //   delete(V);
     //   delete(VU);
	double tmplogDet = log(sqrt(abs(matDet)));
	double tmpFinalRes = -res/2 - tmplogDet;
  //      cout<<"Resulted tmpFinalRes is: "<<exp(-res/2)/sqrt(abs(matDet))<<endl;
	if(tmpFinalRes > 700) 
		return(exp(700));
	return( exp(-res/2)/sqrt(abs(matDet)) );	
}

int PostCal::nextBinary(int * data, int size) {
      //  cout<<"Value of size is: "<<size<<endl;
       // for(long int i=0;i<size;i++)
       //  { 
       //    cout<<" "<<data[i];
       //  }
       // cout<<endl;
	int i = 0;
	int total_one = 0;	
	int index = size-1;
        int one_countinus_in_end = 0;

        while(index >= 0 && data[index] == 1) {
                index = index - 1;
                one_countinus_in_end = one_countinus_in_end + 1;
	}
	if(index >= 0) {
        	while(index >= 0 && data[index] == 0) {
               	 index = index - 1;	
		}
	}
        if(index == -1) {
                while(i <  one_countinus_in_end+1 && i < size) {
                        data[i] = 1;
                        i=i+1;
		}
                i = 0;
                while(i < size-one_countinus_in_end-1) {
                        data[i+one_countinus_in_end+1] = 0;
                        i=i+1;
		}
	}
        else if(one_countinus_in_end == 0) {
                data[index] = 0;
                data[index+1] = 1;
	} else {
                data[index] = 0;
                while(i < one_countinus_in_end + 1) {
                        data[i+index+1] = 1;
			if(i+index+1 >= size)
				printf("ERROR3 %d\n", i+index+1);
                        i=i+1;
		}
                i = 0;
                while(i < size - index - one_countinus_in_end - 2) {
                        data[i+index+one_countinus_in_end+2] = 0;
			if(i+index+one_countinus_in_end+2 >= size) {
				printf("ERROR4 %d\n", i+index+one_countinus_in_end+2);
			}
                        i=i+1;
		}
	}
	i = 0;
	total_one = 0;
	for(i = 0; i < size; i++)
		if(data[i] == 1)
			total_one = total_one + 1;
	
	return(total_one);		
}
int PostCal::decomp(vector<int> &str, int *data, int num)
  {
           int init_start = 0;
 //          cout<<"Come into decomp"<<endl;
           int init_end = 0;
           int test=str.size();
 //          cout<<"Size of the configure is: "<<test<<endl;
           if(str[0]!=-1)
            {
              for (int idx = 0; idx < str.size(); idx++)
               {
                for (int i_x = init_start; i_x<init_start + str[idx]; i_x++)
                 {
                   data[i_x] = 0;
                 }
                int last = init_start + str[idx];
                data[last] = 1;
                init_start = init_start + str[idx] + 1;
              }
             for (int x = init_start + 1; x<num; x++)
              {
                data[x] = 0;
              }
            } else
            {
              for (int y =0; y<num; y++)
               {
                 data[y] = 0;
               }
            }
          int num_x=0;
          for (int x = 0; x<num; x++)
           {
             if(data[x]==1)
              {
                num_x++;
              }
           }
          return num_x;
  }
string PostCal::convert_symbol(int  *data, int num,vector<int> &output)
  {
    string str;
    string sym = "";
    int a = 0;
    int index=0; 
  //  cout<<"Output string infor: "<<endl;
  //  for (int n = 0; n<num; n++)
  //   { 
   //    cout<<" "<<data[n];
   //  }
   // cout<<endl;
    for (int n = 0; n<num; n++)
     {
       if (data[n] == 1)
        {
          ostringstream temp;  //temp as in temporary
          temp<<a;
          string x=temp.str();      //str is temp as string
          output.push_back(a);
          str += x;
          index=1;
          a = 0;
        }
       else
        {
          a++;
        }
     }
    if(index==0)
     {
       output.push_back(-1);
     }
   // for(int i=0; i<output.size(); ++i)
   // cout<< output[i] << ' ';
   // cout<<endl;
   // cout << "The output string is: " << str << endl;
    return str;
  }
double PostCal::computeTotalLikelihood(double * stat, double NCP,map<int,double>& Weight, int nthread) {	
//	int num = 0;
	double sumLikelihood = 0;
	double tmp_likelihood = 0;
	long int total_iteration = 0 ;
	int * configure = (int *) malloc (snpCount * sizeof(int *)); // original data	

	for(long int i = 0; i <= maxCausalSNP; i++)
		total_iteration = total_iteration + nCr(snpCount, i);
	cout << snpCount << endl;
	cout << "Max Causal=" << maxCausalSNP << endl;
	cout << "Total="      << total_iteration << endl;
	for(long int i = 0; i < snpCount; i++) 
		configure[i] = 0;
        //string test[total_iteration];
        vector<vector<int> > test;
       // test.push_back(vector<int>());
        vector<int> res;
        convert_symbol(configure, snpCount,res);
        test.push_back(res);
   //     cout<<"Value of total_iteration is: "<<total_iteration<<endl;
        for(long int i = 1; i < total_iteration; i++) 
         {
        //  test[i]=convert_symbol();
          nextBinary(configure, snpCount);
         // test.push_back(vector<int>());
          vector<int> res;
          convert_symbol(configure, snpCount,res);
          test.push_back(res);
     //     cout<<"It is here "<<i<<endl;
       //   for(long int i=0;i<snpCount;i++)
        //   {
        //     cout<<" "<<configure[i];
        //   }
        //  cout<<endl;
        //  nextBinary(configure, snpCount);          
         }
    //   cout<<"Conversion has finished"<<endl;
       //  cout<<"Hi, I'm here"<<endl;
       // for(long int i = 0; i < total_iteration; i++)
       //  {
       //    for(long int j=0;j<test[i].size();j++)
       //      cout<<" "<<test[i][j];
       //    cout<<endl;
       //  }
       // cout<<endl;
       // delete(configure);
        
       // configure = (int *) malloc (snpCount * sizeof(int *));
        for(long int i=0;i<snpCount;i++)
         {
           configure[i]=0;
         }
        double  prior_x=1;
        for(long int zeng = 0; zeng < snpCount; zeng++)
         {
           double test_x=1;
          // map<int,double>::iterator l_it;
          // l_it = Weight.find(zeng);
          // if(l_it==Weight.end())
          //  {
          //  } else
          //  {
          //    test_x=l_it->second;
          //  }
           if(configure[zeng]==1)
            {
              test_x=test_x;
            } else
            {
              test_x=1-test_x;
            }
           prior_x=prior_x*test_x;
         }
        int test_x=0;
        tmp_likelihood = likelihood(configure, NCP) * prior_x;
        sumLikelihood += tmp_likelihood;
        for(int j = 0; j < snpCount; j++)
         {
           postValues[j] = postValues[j] + tmp_likelihood * configure[j];
         }
      //  histValues[num] = histValues[num] + tmp_likelihood;        
        int confi[nthread][snpCount];
        for(long int i=0;i<nthread;i++)
         {
           for(long int j_x=0;j_x<snpCount;j_x++)
            {
               confi[i][j_x]=configure[j_x];
            }
         }
//        cout<<"Get into parallel region"<<endl;
        int tid;
        omp_set_num_threads(nthread);
        prior_x=1;
        int num;
//        double tmp_likelihood;
   //     cout<<"Get into parallel region"<<endl;
//        #pragma omp parallel for  private (tid,prior_x,num,tmp_likelihood,configure)
        int nloops=0;
   //     cout<<"total_iteration is: "<<total_iteration<<endl;
        #pragma omp parallel   private (tid,prior_x,num,tmp_likelihood,nloops)
          {
             #pragma omp for
	     for(long int i = 1; i < total_iteration; i++) {
   //             cout<<"Already in parallel"<<endl;
                double prior=1;
                tid=omp_get_thread_num();
                num=snpCount;
   //             cout<<"num is: "<<num<<endl;
//                num=decomp(test[i],confi[tid], snpCount);
               #pragma omp critical
               {
               int init_start = 0;
    //           cout<<"Come into decomp"<<endl;
               int init_end = 0;
               
               vector<int> str=test[i];
               
               int test=str.size();
      //         cout<<"Size of the configure is: "<<test<<endl;
               if(str[0]!=-1)
                {
                  for (int idx = 0; idx < str.size(); idx++)
                   {
                     for (int i_x = init_start; i_x<init_start + str[idx]; i_x++)
                       {
                   confi[tid][i_x] = 0;
                 }
                int last = init_start + str[idx];
                confi[tid][last] = 1;
        //        cout<<"value of las is: "<<last;
                init_start = init_start + str[idx] + 1;
              }
        //     cout<<endl;
             for (int x = init_start ; x<num; x++)
              {
                confi[tid][x] = 0;
              }
            } else
            {
              for (int y =0; y<num; y++)
               {
                 confi[tid][y] = 0;
               }
            }
         // int num_x=0;
        //  cout<<"output changed vector"<<endl;
          for (int x = 0; x<num; x++)
           {
            // cout<<"output changed vector"<<endl;
             if(confi[tid][x]==1)
              {
            //    num++;
              //  cout<<confi[tid][x];
              }
 //            cout<<confi[tid][x];
             //cout<<endl;
           }
           }
       
         // return num_x;

          //      cout<<"Value of i is: "<<i<<endl;
                nloops++;
                prior_x=1;
                #pragma omp critical
                 {
                    for(long int z_test = 0; z_test < snpCount; z_test++)
         {
           double test_x=1;
        //   map<int,double>::iterator l_it;
        //   l_it = Weight.find(zeng);
        //   if(l_it==Weight.end())
        //    {
        //    } else
        //    {
        //      test_x=l_it->second;
        //    }
        //   test_x=l_it->second;
           if(confi[tid][z_test]==1)
            {
              test_x=0.01;
            } else
            {
              test_x=1-0.01;
            }
           prior_x=prior_x*test_x;
         }

                 }
          //      cout<<"Value of NCP is: "<<NCP<<endl; 
                tmp_likelihood = likelihood(confi[tid], NCP) * prior_x;
          //      cout<<"tmp_likelihood is: "<<tmp_likelihood<<", prior_x is: "<<prior_x<<endl;
          //      tmp_likelihood = likelihood(configure, NCP) * prior_x;
          //      cout<<"tmp_likelihood is"<<tmp_likelihood<<endl;
                #pragma omp critical
                 {
                   sumLikelihood += tmp_likelihood;
		   for(int j = 0; j < snpCount; j++) 
                    {
                        postValues[j] = postValues[j] + tmp_likelihood * confi[tid][j];
                 //       postValues[j] = postValues[j] + tmp_likelihood * configure[j];
                    }
//		   histValues[num] = histValues[num] + tmp_likelihood;
                 }
 //              cout<<"One test is over"<<endl;
           }
           tid = omp_get_thread_num();

     //      printf("Thread %d performed %d iterations of the loop.\n",tid, nloops );
         }
        free(configure);
   //     cout<<"Parallel is over"<<endl;
        return(sumLikelihood);
}

/*
	stat is the z-scpres
	sigma is the correaltion matrix
	G is the map between snp and the gene (snp, gene)
*/
double PostCal::findOptimalSetGreedy(double * stat, double NCP, char * configure, int *rank,  double inputRho, map<int, double>& Weight, int nthread) {
	int index = 0;
        double rho = 0;
        double total_post = 0;

        totalLikeLihood = computeTotalLikelihood(stat, NCP, Weight, nthread);

	for(int i = 0; i < snpCount; i++)
		total_post += postValues[i];

	printf("Total Likelihood= %e SNP=%d \n", total_post, snpCount);
	
        std::vector<data> items;
        std::set<int>::iterator it;
	//output the poster to files
        for(int i = 0; i < snpCount; i++) {
             //printf("%d==>%e ",i, postValues[i]/total_likelihood);
             items.push_back(data(postValues[i]/total_post, i, 0));
        }
        printf("\n");
        std::sort(items.begin(), items.end(), by_number());
        for(int i = 0; i < snpCount; i++)
                rank[i] = items[i].index1;

        for(int i = 0; i < snpCount; i++)
                configure[i] = '0';
        do{
                rho += postValues[rank[index]]/total_post;
                configure[rank[index]] = '1';
                printf("%d %e\n", rank[index], rho);
                index++;
        } while( rho < inputRho);

        printf("\n");
     //   cout<<"Step of run is over"<<endl;
	return(0);
}
