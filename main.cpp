/*This Code simulates cooling of the Ising model. If you keep the size small
 the output diagram file will allow you to see the lattice. The other file can
 be plotted to see how the magnetization changes in regards to time. In this 
 problem we used the metropolis algorithm.*/

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iomanip>


using namespace std;

double rng(double rngnum);
int row = 15;
int col = 15;
int iterations = 10000;


int main(int argc, char** argv) {
    
    //first we create a lattice of random spins
    double lattice[row][col];
    double time[iterations];
    double magnet[iterations];
    double mag=0;
    
    
    for(int i=0; i<iterations; i++){
        time[i] = i;
    }
    
    double rngind=0;
    double randnum;
    for(int i=0; i<row; i++){
        for(int j=0; j<col; j++){
            randnum = rng(rngind);
            if(randnum < 0.5){
                lattice[i][j] = -1;
            }
            else if(randnum > 0.5){
                lattice[i][j] = 1;
            }
            rngind++;
        }
    }
    
    
    std::remove("Ising_diagrams.txt");
    ofstream diagram;
    std::remove("Ising_data.txt");
    ofstream data;
    
    diagram.open("Ising_diagrams.txt", ios::app);
    for(int i=0; i<row; i++){
        for(int j=0; j<col; j++){
            diagram << (lattice[i][j]+1)/2 << "\t";
        }
        diagram << endl;
    }
    diagram<< endl;
    diagram<< endl;
    diagram.close();
    
    
    
    
      
    //Next we create our energy and location parameters
    double beta = 100;
    double J=3;
    double h=20;
    double E_init, E_flip, E_diff, E_exp;
    double flip_rn;
    int locx, locy;
    
    //Now we implement our metropolis algorithm
    for(int i=0; i<iterations; i++){
        //first we choose a random location
        locx = floor(col*rng(rngind));
        locy = floor(row*rng(rngind+1));
        rngind +=2;
        
        //Now we measure the energy and flipped energy at that location
        //First we consider the boundary conditions
        if(locx == 0 && locy == 0 ){
            E_init = -J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy][locx+1]) - h*lattice[locy][locx];
            E_flip =  J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy][locx+1]) + h*lattice[locy][locx];       
        }
        if(locx == 0 && locy == (row-1) ){
            E_init = -J*(lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx+1]) - h*lattice[locy][locx];
            E_flip =  J*(lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx+1]) + h*lattice[locy][locx];
        }
        if(locx == (col-1) && locy == 0 ){
            E_init = -J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy][locx-1]) - h*lattice[locy][locx];
            E_flip =  J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy][locx-1]) + h*lattice[locy][locx]; 
        }
        if(locx == (col-1)  && locy == (row-1) ){
            E_init = -J*(lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx-1]) - h*lattice[locy][locx];
            E_flip =  J*(lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx-1]) + h*lattice[locy][locx]; 
        }
        if(locx == 0 && locy != 0 && locy != (row-1) ){
            E_init = -J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx+1]) - h*lattice[locy][locx];
            E_flip =  J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx+1]) + h*lattice[locy][locx];  
        }
        if(locx == (col-1)  && locy != 0  && locy != (row-1) ){
            E_init = -J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx-1]) - h*lattice[locy][locx];
            E_flip =  J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx-1]) + h*lattice[locy][locx];  
        }
        if(locx != 0 && locx != (col-1)  && locy == 0 ){
            E_init = -J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy][locx+1] + lattice[locy][locx]*lattice[locy][locx-1]) - h*lattice[locy][locx];
            E_flip =  J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy][locx+1] + lattice[locy][locx]*lattice[locy][locx-1]) + h*lattice[locy][locx]; 
        }
        if(locx != 0 && locx != (col-1)  && locy == (row-1) ){
            E_init = -J*(lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx+1] + lattice[locy][locx]*lattice[locy][locx-1]) - h*lattice[locy][locx];
            E_flip =  J*(lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx+1] + lattice[locy][locx]*lattice[locy][locx-1]) + h*lattice[locy][locx]; 
        }
        //And then the general case
        else {
            E_init = -J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx+1] + lattice[locy][locx]*lattice[locy][locx-1]) - h*lattice[locy][locx];
            E_flip =  J*(lattice[locy][locx]*lattice[locy+1][locx] + lattice[locy][locx]*lattice[locy-1][locx] + lattice[locy][locx]*lattice[locy][locx+1] + lattice[locy][locx]*lattice[locy][locx-1]) + h*lattice[locy][locx];
        }
        
        //We then calculate the difference in energy
        E_diff = E_flip - E_init;
        
        //We get another random number
        flip_rn = rng(rngind);
        rngind++;
        //And we calculate our exponential energy
        E_exp = exp(-beta*E_diff);
        //We then see if we flip or not
        if(flip_rn <= E_exp ){
            lattice[locy][locx] *= -1;
        }
        mag=0;
        for(int i=0; i<row; i++){
            for(int j=0; j<col; j++){
                mag += lattice[i][j];
            }
        }
        
        data.open("Ising_data.txt", ios::app);
        data << time[i] << "\t" << mag << endl;
        data.close();
    }
    
    
    diagram.open("Ising_diagrams.txt", ios::app);
    for(int i=0; i<row; i++){
        for(int j=0; j<col; j++){
            diagram << (lattice[i][j]+1)/2 << "\t";
        }
        diagram << endl;
    }
    diagram<< endl;
    diagram<< endl;
    diagram.close();
    
    
    
    return 0;
}

double rng (double rngnum){
    srand(time(0) + rngnum);
    double epsilon = (double) rand()/RAND_MAX;
    return epsilon;
}