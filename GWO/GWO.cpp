#include "GWO.hpp"

using namespace std;
int len(int x[]) 
{
	return (sizeof(x)/sizeof(x[0]));
}
int lenDouble(double x[])
{
	return (sizeof(x)/sizeof(x[0]));
}
double *power(double x[], int y) 
{
	int lx=lenDouble(x);
	double *xp= new double[lx]; 
	for(int i=0;i<lx;i++)
	{
		xp[i]=pow(x[i],y);
	}
	return (xp);
}
double SummedUfun(double x[],int a,int k,int m)
{
	int dim=lenDouble(x);
	double y[dim];
	for (int i=0;i<dim;i++)
	{	
		if(x[i]>a)
		{
			y[i]=k*pow((x[i]-a),m);
		}
		else if(x[i]<-a)
		{
			y[i]=k*pow((-x[i]-a),m);
		}
		else
		{
			y[i]=0;
		}
	}
	double summedy=0;
	for (int i=0;i<dim;i++)
	{
		summedy=summedy+y[i];
	}
	return (summedy);
}
double prod(double a[])
{
	int l=lenDouble(a);
	double result;
	for (int i=0;i<l;i++)
	{
		result=a[i]*result;
	}
	return result;
}
double arraySum(double a[])
{
	int n = lenDouble(a); 
    double sum  = 0;
    for(int i=0;i<n;i++)
    {
    	sum=sum+a[i];
    }  
    return (sum); 
}
double F1(double x[])
{
	double fit=arraySum(power(x,2));
	return fit;
}
double F2(double x[])
{	
	double fit=arraySum(x)+prod(x);
	return fit;
}
double max(double x[])
{
	int l=lenDouble(x);
	return ( *std::max_element(x, x+l) );
}
double F3(double x[])
{
	double dim=lenDouble(x)+1;
	double fit=0;
	int i;
	for(i=1;i<dim;i++)
	{
		double sliced[i];
		for (int j=0; j<i+1;j++)
		{
			sliced[j]=x[j];
		}
		fit=fit+pow(arraySum(sliced),2);
	}
	return fit;
}
double F4(double x[])
{
	double fit=abs(max(x));
	return fit;
}
double F5(double x[])
{
	int dim=lenDouble(x);
	double sliced1[dim];
	double sliced2[dim];
	for (int i=1;i<dim+1;i++)
	{	
		sliced1[i]=x[i];
		sliced2[i-1]=pow(x[i-1],2);
	}
	double slicedpower[dim];
	for (int i=0;i<dim;i++)
	{
		slicedpower[i]=100*(sliced1[i]-sliced2[i]);
	}
	double slicedpowerpower[dim];
	for (int i=0;i<dim;i++)
	{
		slicedpowerpower[i]=pow(slicedpower[i],2);
	}
	double sliced3powered[dim];
	for(int i=0;i<dim;i++)
	{
		sliced3powered[i]=pow(x[i]-1,2);
	}
	double finalarray[dim];
	for (int i=0;i<dim;i++)
	{
		finalarray[i]=slicedpowerpower[i]+sliced3powered[i];
	}
	double fit=arraySum(finalarray);  
	return fit;
}
double F6(double x[])
{	
	int dim=lenDouble(x);
	double chngx[dim];
	for(int i=0;i<dim;i++) 
	{
		chngx[i]=x[i]+0.5;
	}
	double fit=pow(arraySum(chngx),2);
	return fit;
}
double F7(double x[])
{
	double fit=0;
	int dim=lenDouble(x);
	for(int i=0;i<dim;i++)
	{
		fit=fit+((i+1)*pow(x[i],4)+rand()%2); 
	}
	return fit;
}
double F8(double x[])
{
	int dim=lenDouble(x);
	double arraypowered[dim];
	for (int i=0;i<dim;i++)
	{
		arraypowered[i]=x[i]*sin(pow(x[i],1/2));
	}
	double fit=-(arraySum(arraypowered));
	return fit;
}
double F9(double x[])
{
	int dim=lenDouble(x);
	double arraypowered[dim];
	for (int i=0;i<dim;i++)
	{
		arraypowered[i]=pow(x[i],2)-10*cos(2*M_PI*x[i]);
	}
	double fit=arraySum(arraypowered)+10*dim;
	return fit;
}
double F10(double x[])
{	
	int dim=lenDouble(x);
	double *arraypowered=power(x,2);
	double arraysumrooted=pow(arraySum(arraypowered)/dim,1/2);
	double arraycosed[dim];
	for (int i=0;i<dim;i++)
	{
		arraycosed[i]=cos(2*M_PI*x[i]);
	} 
	double arraysumofcosed=arraySum(arraycosed);
	double fit=-20*exp(-0.2*arraysumrooted)-exp(arraysumofcosed/dim)+20+exp(1);
	return fit;
}
double F11(double x[]) 
{
	int dim=lenDouble(x);
	double *arraysquared=power(x,2); 
	double summation=arraySum(arraysquared);
	double production=1;
	for (int i=1;i<dim;i++)
	{
		production=production*(cos(x[i-1]/pow(i,1/2)));
	}
	float fit=(summation)/4000-production+1;
	return fit;
}
double F12(double x[])
{
	int dim=lenDouble(x);
	double slicedAndDividedbyFour[dim];
	for (int i=1;i<dim-1;i++)
	{
		slicedAndDividedbyFour[i]=(x[i]+1)/4;
	}
	double slicedAndDividedbyFouraddOne[dim];
	for (int i=0;i<dim;i++)
	{
		slicedAndDividedbyFouraddOne[i]=1+10*pow(sin(M_PI*(1+slicedAndDividedbyFour[i])),2);
	}
	double finalSummation;
	for (int i=0;i<dim;i++)
	{
		finalSummation=finalSummation+(pow(slicedAndDividedbyFour[i],2)*
								slicedAndDividedbyFouraddOne[i]);
	}
	double fit=(M_PI/dim)*(10*pow((sin(M_PI*(1+(x[0]+1/4)))),
		2)+finalSummation+pow(((x[dim-1]+1)/4),
		2))+SummedUfun(x,10,100,4);
	return fit;
}
double F13(double x[])
{
	int dim=lenDouble(x);
	double slicedminusone[dim];
	for(int i=0;i<dim-1;i++)
	{
		slicedminusone[i]=pow(x[i]-1,2);
	}
	double slicedfrom1[dim];
	for (int i=0;i<dim;i++)
	{
		slicedfrom1[i]=slicedminusone[i]*(1+pow(sin(3*M_PI*x[i]),2));
	}
	double fit=0.1*(pow(sin(3*M_PI*x[1]),2)+arraySum(slicedfrom1)+
		pow(x[dim-1]-1,2)*(1+pow(sin(2*M_PI*x[dim-1]),2)))+SummedUfun(x,5,100,4);
	return fit;
}
double F14(double x[])
{
	double aS[][25]={{-32,-16,0,16,32,-32,-16,0,16,32,-32,-16,0,
		16,32,-32,-16,0,16,32,-32,-16,0,16,32},{-32,-32,-32,
		-32,-32,-16,-16,-16,-16,-16,0,0,0,0,0,16,16,16,16,16,
		32,32,32,32,32},};
	double bS[25],H[2];
	for(int i=0;i<25;i++)
	{
		for (int j=0;j<2;j++)
		{
			H[j]=x[i]-aS[j][i];
		}
		bS[i]=arraySum(power(H,6));
	}
	double wplusBs[25];
	for (int i=0;i<24;i++)
	{
		wplusBs[i]=1/(i+1+bS[i]);
	}
	wplusBs[25]=25;
	double fit=(1/500)+pow(arraySum(wplusBs),-1);
	return (fit);
}

int GWO(int lb,int ub,int dim,int SearchAgents_no,int Max_iter,std::ofstream& myFile)
{
	double Alpha_pos[dim];
	double Alpha_score = std::numeric_limits<double>::infinity();
	double Beta_pos[dim];
	double Beta_score = std::numeric_limits<double>::infinity();
	double Delta_pos[dim];
	double Delta_score = std::numeric_limits<double>::infinity();
	double Positions[SearchAgents_no][dim];
    std::random_device rd;
    std::mt19937_64 mt(rd());
    std::uniform_real_distribution<double> distribution(lb, ub);
    double array[SearchAgents_no][dim];
    for(int i = 0; i < SearchAgents_no; i++)
    {
    	for (int j=0;j<dim;j++)
    	{
        	double d = distribution(mt);
        	array[i][j] = d;
        }
    }
	double productArray[SearchAgents_no][dim];
	for (int i=0;i<SearchAgents_no;i++)
	{
		for(int j=0;j<dim;j++)
		{
			productArray[i][j]=array[i][j]*(ub-lb);
		}
	}
	double addArray[SearchAgents_no][dim];
	for (int i=0; i<SearchAgents_no;i++)
	{
		for (int j=0;j<dim;j++)
		{
			addArray[i][j]=productArray[i][j]+lb;
		}
	}
	for (int i=0;i<SearchAgents_no;i++) 
	{																			
		for(int j=0;j<dim;j++)
		{
			Positions[i][j]=addArray[i][j];												
		}
	}
	double Convergence_curve[Max_iter];
	double fitness;
	for(int l=0;l<Max_iter;l++)
	{
		for(int i=0;i<SearchAgents_no;i++)
		{	
			double fitness=F14(Positions[i]);	
			if(fitness<Alpha_score)
			{
				Alpha_score=fitness;
				for(int z=0;z<SearchAgents_no;z++)
				{
					Alpha_pos[z]=Positions[i][z];
				}
			}
			if(fitness>Alpha_score && fitness<Beta_score)
			{
				Beta_score=fitness;
				for(int z=0;z<SearchAgents_no;z++)
				{
					Beta_pos[z]=Positions[i][z];
				}
			}
			if(fitness>Alpha_score && fitness>Beta_score && fitness<Delta_score)
			{
				Delta_score=fitness;
				for(int z=0;z<SearchAgents_no;z++)
				{
					Delta_pos[z]=Positions[i][z];
				}
			}
		}
		int a=2-1*((2)/Max_iter);
		for (int i=0;i<SearchAgents_no;i++)
		{
			for (int j=0;j<dim;j++)
			{
				int r1=rand()%2;
				int r2=rand()%2;
				int A1=2*a*r1-a;
				int C1=2*r2;
				int D_alpha=abs(C1*Alpha_pos[j]-Positions[i][j]);
				int X1=Alpha_pos[j]-A1*D_alpha;

				r1=rand()%2;
				r2=rand()%2;
				int A2=2*a*r1-a;
				int C2=2*r2;
				int D_beta=abs(C2*Beta_pos[j]-Positions[i][j]);
				int X2=Beta_pos[j]-A2*D_beta;
				r1=rand()%2;
				r2=rand()%2;
				int A3=2*a*r1-a;
				int C3=2*r2;
				int D_delta=abs(C2*Delta_pos[j]-Positions[i][j]);
				int X3=Delta_pos[j]-A3*D_delta;
				Positions[i][j]=(X1+X2+X3)/3;
			}
		}
		myFile <<Alpha_score<<',';
	}
	myFile<<endl;
}