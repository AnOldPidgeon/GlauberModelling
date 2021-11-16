// This program should make graphs that replicate Fig. 4
// in this paper about the Glauber Model:
// http://nuclear.ucdavis.edu/~calderon/nuclearDocs/groupMeeting/glauberModelAnnReview09052007.pdf
// but model Pb<->Pb collision instead of Au<->Au

#include "math.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TF1.h"
#include "TRandom3.h"
#include "TVector3.h"
#include "TRotation.h"           // I probably went overboard with the things I've imported here
#include "TLorentzVector.h"      // BUT YOU CAN'T BE TOO SURE
#include <iostream>
#include <cmath>

//here, we are plotting out lead nucleons of colliing atoms on the ZX and XY Plane

void collisions2 () {      //rename this at your convenience
    
    
    // Woods-Saxon distribution variables
    
    int a_nuc = 208;   // number of nucleons in one nucleus, should be 208 for lead
    double a_sd = 0.54;   // skin depth
    double r_max = 6.34;  // used to calculate a_sd, not really used in the code
    
    int events = 100; // the number of collisions of lead atoms, also size of
    
    //boolean values. I hate these:(
    int loopAPart = 0; //boolean value. return one under right condition
    int loopBPart = 0;
    
    //collision indices. to be added to respective boolean values ugh
    
    int AIndex = 0;
    int BIndex = 0;
    
    //for counting participating nucleons
    int APart = 0; // number of participating nucleons in part A
    int BPart = 0;
    
    double pi = 3.14159265358979323; // yes, this precision is important lol
    TRandom3 rnd3(pi);     // initialise random number generator with unique seed
    
    double d = sqrt(42.0/(10.0*pi));// 42 mb at sqrt(s)=200 GeV,  64 mb at sqrt(s)=2.76 TeV, and 70 mb at sqrt(s)=5.02 TeV
    
    
    // this will give us the spherical distribution of nucleons,
    //  based on the Woods-Saxon distribution equation.
    TF1* get_rho = new TF1("get_rho", "x*x*(1) / (1 + exp((x-6.34)/0.54))", 0, 10); // to generate the random values for radius
    TF1* get_theta = new TF1("get_theta", "sin(x)", 0, pi);
    
    //make arrays of booleans
    
    int ABool[208] = {};
    int BBool[208] = {};
    
    // Make the arrays for storing nucleon A and B positions
    double x_A[208] = {};
    double y_A[208] = {};
    double z_A[208] = {}; // x, y, and z positions for all nuclons in nucleus A
    
    double x_AColl[208] = {}; // we only care about the coordinates in x and y space for colliding nucleons???
    double y_AColl[208] = {};
    double z_AColl[208] = {};
    
    double x_B[208] = {};
    double y_B[208] = {};
    double z_B[208] = {}; // x, y, and z positions for all nuclons in nucleus B
    
    double x_BColl[208] = {};
    double y_BColl[208] = {};
    double z_BColl[208] = {};
    
    int nColl = 0; //total collisions of one nucleon
    int eventColl = 0; // total nuber of collisions in the event
    
    double collArr[100] = {}; // to store the participating nucleons in each collision. Size determined by num of events
    
    
    // Make variables for spherical coordinates
    double r;
    double theta;
    double phi;
    
    
    
    // loop for finding nucleus A's nucleon positions
    for (int i = 0; i < a_nuc; i++)
    {
        // randomly generate a radius away from the centre of the nucleus, according to the Woods-Saxon distribution
        r = get_rho->GetRandom();
        theta = get_theta->GetRandom();       // separately generate a random angle theta based on sine distribution, from 0 to pi
        phi = 2*pi*rnd3.Rndm();               // separately generate an angle phi, uniformly in the range from 0 to 2pi
        
        x_A[i] = r * cos(phi) * sin(theta) - 3;   // find the x,y,z coordinates, convert from spherical variables to cartesian, and subtract by half the inpact parameter.
        y_A[i] = r * sin(phi) * sin(theta);
        z_A[i] = r * cos(theta)-12; // z axis is the beam axis
        
    } // end of loop
    
    
    
    // loop for finding nucleus B's nucleon positions
    for (int i = 0; i < a_nuc; i++)
    {
        r = get_rho->GetRandom();
        theta = get_theta->GetRandom();
        phi = 2*pi*rnd3.Rndm();
        
        x_B[i] = r * cos(phi) * sin(theta) + 3; // Shift the x coordinate to the right by 3 fermi (half the impact parameter)
        y_B[i] = r * sin(phi) * sin(theta);
        z_B[i] = r * cos(theta)+12;
        
    } // end of loop
    
    // Count the number of participating nucleons in the collision (for now only one iteration)
    
    for (int j = 0; j < a_nuc; j++){ //for each nucleon in nucleus A
        
        for(int k = 0; k < a_nuc; k++){ //counting nucleons in nucleus B
            
            if(sqrt(pow(x_A[j]-x_B[k], 2)+pow(y_A[j]-y_B[k], 2)) <= d){ // continue loop iff this condition is met
                
                nColl++; //check for each iteration that both radii are less than or equal to the distance
                
                x_AColl[AIndex] = x_A[j]; //assigning values to collision array
                y_AColl[AIndex] = y_A[j];
                z_AColl[AIndex] = z_A[j];
                //include z coordinate!
                
                x_BColl[BIndex] = x_B[k];
                y_BColl[BIndex] = y_B[k];
                z_BColl[BIndex] = z_B[k];
                //here too
                
                ABool[j] = 1; //return 0 otherwise boi
                BBool[k] = 1;
            }
            
            /*if(sqrt(pow(x_A[k]-x_B[j], 2)+pow(y_A[k]-y_B[j], 2)) <= d){ // just switched a bunch of indices
                
                x_BColl[BIndex] = x_B[j];
                y_BColl[BIndex] = y_B[j];
                //here too
                
                loopBPart = 1;
            }*/
        }
        
        AIndex += ABool[j]; // only increments if the above if statement is true, should be 0 otherwise
        BIndex += ABool[j];
    }
    // end of loop
    
    for (int n = 0; n < 208; n++){ //to count the number of participating nucleons
        
        APart += ABool[n];
        BPart += BBool[n];
    }
        
    cout << "The number of participating nucleons is: " << APart + BPart << endl;
    cout << "The total number of collisions is: " << nColl << endl;
    
    // MAKING THE 3D GRAPH, I think I fucked up here
    
    /*TCanvas* c1 = new TCanvas ("c1", "3D Pb Nucleus A", 500, 500, 500, 500);
    
    TGraph2D* nucA = new TGraph2D(a_nuc, x_A, y_A, z_A); // this is a 3D graph of nucleus A
    nucA-> SetTitle("3D Distributions; x; y; z");
    nucA-> SetMarkerStyle(8);
    nucA-> SetMarkerSize(2);
    nucA-> SetMarkerColor(kOrange+6);
    nucA-> Draw("AP");
    
    TGraph2D* nucAColl = new TGraph2D(a_nuc, x_AColl, y_AColl, z_A); // this is a 3D graph of colliding nucleons in nucleus A
    nucAColl-> SetMarkerStyle(8);
    nucAColl-> SetMarkerSize(2);
    nucAColl-> SetMarkerColor(kRed);
    nucAColl-> Draw("P SAME");
    
    //TCanvas* c2 = new TCanvas ("c2", "3D Pb Nucleus B", 500, 500, 500, 500);
    TGraph2D* nucB = new TGraph2D(a_nuc, x_B, y_B, z_B); // this is a 3D graph of nucleus B
    nucB-> SetMarkerStyle(8);
    nucB-> SetMarkerSize(2);
    nucB-> SetMarkerColor(kAzure+6);
    nucB-> Draw("P SAME");  // this means it will be on the same graph as nucleus A
    
    TGraph2D* nucBColl = new TGraph2D(a_nuc, x_BColl, y_BColl, z_B); // this is a 3D graph of colliding nucleons in nucleus B
    nucBColl-> SetMarkerStyle(8);
    nucBColl-> SetMarkerSize(2);
    nucBColl-> SetMarkerColor(kBlue);
    nucBColl-> Draw("P SAME");*/
    
    // MAKING 2D GRAPH IN XZ PLANE
    
    TCanvas* c2 = new TCanvas ("c2", "2D Pb Nucleus A", 500, 500, 500, 500);
    
    TH2D* modelXZ = new TH2D("histo2d","2D Distributions in ZX space; z (fm); x (fm)",2000,-20,20,2000,-10,10); //Dummy histogram
    modelXZ->Draw();
    
    /*TArrow *ar3 = new TArrow(-12,-3,0,-3,1,"|>");
    ar3->SetAngle(90);
    ar3->SetLineWidth(2);
    ar3->Draw("P SAME");*/
    
    TGraph* xz_nucA = new TGraph(a_nuc, z_A, x_A);  // this is a 2D graph of nucleus A
    //xz_nucA->GetXaxis()->SetRangeUser(-24,24);
    //xz_nucA->SetTitle("2D Distributions in XZ space; z (fm); x (fm)");
    xz_nucA->SetMarkerStyle(8);
    xz_nucA->SetMarkerSize(2);
    xz_nucA-> SetMarkerColor(kOrange+6);
    xz_nucA->Draw("P SAME");
    
    TGraph* xz_nucAColl = new TGraph(APart, z_A, x_AColl);  // this is a 2D graph of colliding nucleons in nucleus A
    xz_nucAColl->SetMarkerStyle(8);
    xz_nucAColl->SetMarkerSize(2);
    xz_nucAColl->SetMarkerColor(kRed);
    xz_nucAColl->Draw("P SAME");
    
    TGraph* xz_nucB = new TGraph(a_nuc, z_B, x_B);  // this is a 2D graph of nucleus B
    xz_nucB->SetMarkerStyle(8);
    xz_nucB->SetMarkerSize(2);
    xz_nucB->SetMarkerColor(kAzure+6);
    xz_nucB->Draw("P SAME");  // this means it will be on the same graph as nucleus A
    
    TGraph* xz_nucBColl = new TGraph(BPart, z_B, x_BColl);  // this is a 2D graph of colliding nucleons in nucleus B
    xz_nucBColl->SetMarkerStyle(8);
    xz_nucBColl->SetMarkerSize(2);
    xz_nucBColl->SetMarkerColor(kBlue);
    xz_nucBColl->Draw("P SAME");
    
    // MAKING THE 2D GRAPH IN XY PLANE
    
    TCanvas* c3 = new TCanvas ("c3", "2D Pb Nucleus A", 500, 500, 500, 500);
    
    TH2D* modelXY = new TH2D("histo2d","2D Distributions in XY xpace; x (fm); y(fm)",2000,-15,15,2000,-10,10); //Dummy histogram
    modelXY->Draw();
    
    TGraph* flat_nucA = new TGraph(a_nuc, x_A, y_A);  // this is a 2D graph of nucleus A
    //flat_nucA->SetTitle("2D Distributions in XY xpace; x (fm); y (fm)");
    flat_nucA->SetMarkerStyle(8);
    flat_nucA->SetMarkerSize(2);
    flat_nucA-> SetMarkerColor(kOrange+6);
    flat_nucA->Draw("P SAME");
    
    TGraph* flat_nucAColl = new TGraph(a_nuc, x_AColl, y_AColl);  // this is a 2D graph of colliding nucleons in nucleus A
    flat_nucAColl->SetMarkerStyle(8);
    flat_nucAColl->SetMarkerSize(2);
    flat_nucAColl->SetMarkerColor(kRed);
    flat_nucAColl->Draw("P SAME");
    
    TGraph* flat_nucB = new TGraph(a_nuc, x_B, y_B);  // this is a 2D graph of nucleus B
    flat_nucB->SetMarkerStyle(8);
    flat_nucB->SetMarkerSize(2);
    flat_nucB->SetMarkerColor(kAzure+6);
    flat_nucB->Draw("P SAME");  // this means it will be on the same graph as nucleus A
    
    TGraph* flat_nucBColl = new TGraph(a_nuc, x_BColl, y_BColl);  // this is a 2D graph of colliding nucleons in nucleus B
    flat_nucBColl->SetMarkerStyle(8);
    flat_nucBColl->SetMarkerSize(2);
    flat_nucBColl->SetMarkerColor(kBlue);
    flat_nucBColl->Draw("P SAME");
    
    
} // end of program#include "collisions2.h"
