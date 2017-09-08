// To compile: g++ CosSq.cpp -o CosSq.out

#include <vector>
#include <iostream>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <time.h> 
using namespace std;

class Atom{
public:
double x,y,z;
};

int main(int argc, char* argv[])
{
	ifstream input;
	ofstream log;
	ofstream output; //<cos>, <cos^2>, ...
	ofstream output2; //Dipole polarization of H2O in z direction
	input.open(argv[1]);
	output.open(argv[2]);
	output2.open(argv[3]);
	log.open(argv[4]);

	time_t nStart = time(NULL);//Start timer
	
	int nmol=2654; // # of water molecules in gro file
	int nbin=800; // # of slices in z direction
	int skip=0; // # of ps to be skipped (1ps might contain 10 frames, check mdp settings)
	double P0=0.048946478141643074; //Dipole moment of a SPC/E H2O molecule
	int dbg=0;

	string line; // to store each line scanned from the input
	string line2[3];// to store ow, hw1, and hw2 strings
	vector<Atom> ow_ary, hw1_ary, hw2_ary;// collect x,y,z of atoms from line2; vector size = nmol; clear every frame
	double box_x, box_y, box_z, dz, dV;

	//Obtain box dimensions
	for ( int i=0; i<nmol*3+3; i++ )
	{
		getline(input,line);
	}
	box_x = atof(line.substr(0,10).c_str());
	box_y = atof(line.substr(11,10).c_str());
	box_z = atof(line.substr(21,10).c_str());
	dz = box_z/nbin;
	dV= box_x*box_y*dz;
	input.clear();
	input.seekg(0,ios::beg); //reset flag
	
	//Declare an array to store number of molecules in each bin in each frame
	int *nmolbin;
	nmolbin=(int *)malloc(nbin * sizeof(int)); //array to store # of atoms in each bin

        //Declare arrays to store sum of cos, cos^2, cos^3, and cos^4 in each bin (Zeros every time frame)
        double *cosbin, *cos2bin, *cos3bin, *cos4bin, *dippol;
        cosbin=(double *)malloc(nbin * sizeof(double));
        cos2bin=(double *)malloc(nbin * sizeof(double));
        cos3bin=(double *)malloc(nbin * sizeof(double));
        cos4bin=(double *)malloc(nbin * sizeof(double));
	dippol=(double *)malloc(nbin * sizeof(double));

        for (int i=0; i<nbin; i++)
        {
        	nmolbin[i]=0;
		dippol[i]=0;
		cosbin[i]=0; cos2bin[i]=0; cos3bin[i]=0; cos4bin[i]=0;
        }


	//Skipping unuseful frames, i.e. first 500ps   <<<<<<Subroutine>>>>>>
	int nskipframe=0;
	for (int s=0; s<skip*10; s++)
	{
		for (int ss=0; ss<nmol*3+3; ss++)
		{
			getline(input,line);
		}
		
		nskipframe++;
		if (nskipframe % 200 ==0)
		{
			cout << "\rSkipping frame " << nskipframe << "       ";
			cout << "Time(ps) " << nskipframe/10 << flush;
			log << "Skipping frame " << nskipframe << "       ";
			log << "Time(ps) " << nskipframe/10 << endl;
		}
		if (s == skip*10-1)
		{
			cout << endl;
		}
	}

	//Starting scanning molecular orientations
	int nframe=0;//Initialize counting # of frame
	int nbin_thread;
	double norm[3] = {0,0,1}, dipole[3], cosine;
	bool finish_flag=false;
	#pragma omp parallel num_threads(4) private(nbin_thread,dipole,cosine)
	{	
		//while(!input.eof())
		while(finish_flag==false)
		{
			#pragma omp single
			{	
				getline(input,line); //skip 1st line: time inform
				if (line.empty())
				{
					//break; //breaking out loop forbidden in OpenMP block
					finish_flag=true;
					cout << "pass" << endl;
				}
				else
				{
					nframe++; //count # of frames
					//Print reading progress to screen (every 100 frames)
					if (nframe % 200 == 0)
					{
						cout << "\rReading frame " << nframe << "      ";
						cout << "Time(ps) " << nframe/10 << flush;
						log << "Reading frame " << nframe << "      ";
						log << "Time(ps) " << nframe/10 << endl;
					}

					getline(input,line); //skip 2nd line: # of atoms
				}

				if(finish_flag==false)
				{
					Atom atom;
					for (int m=0; m<nmol; m++) // For every molecule in this frame
					{
						for (int n=0; n<3; n++) // For evey atom in this molecule. H2O contains 3 atoms
						{
							getline(input,line2[n]);
						}
						atom.x = atof(line2[0].substr(20,8).c_str());
						atom.y = atof(line2[0].substr(28,8).c_str());
						atom.z = atof(line2[0].substr(36,8).c_str());
						ow_ary.push_back(atom);
						atom.x = atof(line2[1].substr(20,8).c_str());
						atom.y = atof(line2[1].substr(28,8).c_str());
						atom.z = atof(line2[1].substr(36,8).c_str());
						hw1_ary.push_back(atom);
						atom.x = atof(line2[2].substr(20,8).c_str());
						atom.y = atof(line2[2].substr(28,8).c_str());
						atom.z = atof(line2[2].substr(36,8).c_str());
						hw2_ary.push_back(atom);	
					}

					//Make sure data are read correctly   <<<<<<Subroutine>>>>>>
					if ( (ow_ary.size() != nmol ) || (hw1_ary.size() !=nmol) || (hw2_ary.size() != nmol) )
					{
						cout << "Atom numbers error! End Program." << endl;
						cout << "Input molecule number: " << nmol << endl;
						cout << "Size of ow_ary: " << ow_ary.size() << endl;
						cout << "Size of hw1_ary: " << hw1_ary.size() << endl;
						cout << "Size of hw2_ary: " << hw2_ary.size() << endl;
						//exit (EXIT_FAILURE);
						finish_flag=true;
					}
				}
			}	

			//Initialize the sum of molecule numbers, cos, cos^2, etc. in each bin with value 0. The arrays will be zero every frame.
			if(finish_flag==false)
			{
				#pragma omp for
				for (int j=0; j<nmol; j++)
				{
					dipole[0] = (hw1_ary[j].x+hw2_ary[j].x)/2 - ow_ary[j].x;
					dipole[1] = (hw1_ary[j].y+hw2_ary[j].y)/2 - ow_ary[j].y;
					dipole[2] = (hw1_ary[j].z+hw2_ary[j].z)/2 - ow_ary[j].z;
					cosine = ( norm[0]*dipole[0]+norm[1]*dipole[1]+norm[2]*dipole[2] ) / ( sqrt(dipole[0]*dipole[0]+dipole[1]*dipole[1]+dipole[2]*dipole[2])*(1) );
					nbin_thread = min(floor(ow_ary[j].z/dz),(double)nbin);

					#pragma omp critical
					{
						nmolbin[nbin_thread]++;
						cosbin[nbin_thread] += cosine;
						cos2bin[nbin_thread] += cosine*cosine;
						cos3bin[nbin_thread] += cosine*cosine*cosine;
						cos4bin[nbin_thread] += cosine*cosine*cosine*cosine;
					}
				}
				#pragma omp barrier //not required in this case

				#pragma omp single
				{	
					ow_ary.clear(); //empty atomic vectors before starting scanning the next frame
					hw1_ary.clear();
					hw2_ary.clear();
					getline(input,line); //skip the last line of each frame: box dimensions
				}
			}
		}//<--End of the loop over time frames
	}//<-- omp forking finish 


	output << "@	title: Water orientation with respect to normal" << endl;
        output << "@    xaxis: box (nm)" << endl;
        output << "@    yaxis: <cosine>, <cosine^2>, <cosine^3>, <cosine^4>" << endl;
        output2 << "@    title: Water dipole polarization" << endl;
        output2 << "@    xaxis: box (nm)" << endl;
        output2 << "@    yaxis: P_dip (e/nm^2)" << endl;

	//Devide the sum of the number-averaged cos, cos^2, .. by total number of time frames
	for (int i=0; i<nbin; i++)
	{
                dippol[i] = cosbin[i]*P0/(dV)/( (double) nframe);
                if (nmolbin[i] !=0)
                {
                        cosbin[i] = cosbin[i]/nmolbin[i]; cos2bin[i] = cos2bin[i]/nmolbin[i];
                        cos3bin[i] = cos3bin[i]/nmolbin[i]; cos4bin[i] = cos4bin[i]/nmolbin[i];
                }
                else
                {
                        cosbin[i]=0; cos2bin[i]=0; cos3bin[i]=0; cos4bin[i]=0;
                }

		output << setw(18) << setprecision(6) << fixed << box_z/nbin*i;
		output << setw(18) << setprecision(6) << fixed << cosbin[i]; 
                output << setw(18) << setprecision(6) << fixed << cos2bin[i];
                output << setw(18) << setprecision(6) << fixed << cos3bin[i];
                output << setw(18) << setprecision(6) << fixed << cos4bin[i] << endl;
		output2 << setw(18) << setprecision(6) << fixed << box_z/nbin*i;
		output2 << setw(18) << setprecision(6) << fixed << dippol[i] << endl;
	}	

	//Estimate cpu time for the calculations
	time_t nEnd = time(NULL);//End timer
	cout<<"Elapsed time is :  "<< nEnd-nStart << " seconds " << endl;
	log << "Elapsed time is :  "<< nEnd-nStart << " seconds " << endl;

	//Output debug information
	if (dbg==0)
	{
		int sum=0;
		log << "\nThe bgin of debug information:" << endl;
		for (int i=0; i<nbin; i++)
		{
			log << "Number of SOL molecules in " << i << "th bin: "  << nmolbin[i] << endl;
			sum += nmolbin[i];//Total number of molecules (sum of frames)
		}
		log << "\nTotal SOL molecules: " << sum << endl;
		log << "The end of begug information" << endl;
	}

	//Relase memory for dynamic arrays
	free(nmolbin); free(cosbin); free(cos2bin); free(cos3bin); free(cos4bin); free(dippol);

	//Close files
	input.close();
	output.close();
	output2.close();
	log.close();

	return 0;
}
