// Symplectic.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "Symplectic.hpp"

int main(int argc, char** argv)
{
    /*============================================================================================
     BEGIN PULSE INITIALISATION
    ============================================================================================*/
    //turns out the largest size you can do at a time is r=12
    const unsigned int rBATCH = 12;
    const unsigned int rTOTAL = 20;
    const unsigned int sBATCH = pow(2, rBATCH); //number of particles per batch (size of batch)
    const unsigned int sTOTAL = pow(2, rTOTAL); //total size
    const unsigned int nBATCH = pow(2, rTOTAL - rBATCH); //number of batches
    Set beam(sTOTAL);
    beam.SetZdist(1.60, 3.1398e-5, 3, 3);
    beam.SetRXYPhidist(0, 4.972e-4, 3);
    beam.SetGBRXYPhidist(1.0074, 8e-4, 4);
    const double G = 1.2344; 
    const double GB = std::sqrt(G*G - 1);
    beam.SetGBZdist(GB, 0); //check definition tomake sure no div by 0
    beam.SetScale(1.2, 1.0, 1.0, 1.0, 1.0, 1.0);
    beam.AddZdiv(1.60, 1e5);
    beam.AddXdiv(0, 15);
    beam.AddYdiv(0, 15);
    beam.SetGBXEmittance(1.2 * 1.2398e-7);
    beam.SetGBYEmittance(1.2 * 1.2398e-7);
    Set beamOUT = beam;
    //write_current_to_file("C:/Users/SimulationWorkstatio/Desktop/Symplectic", beam, 0);

    /*============================================================================================
    END INITIALISATION
    ============================================================================================*/

    /*============================================================================================
    BEGIN MATRIX MAPPING INITIALISATION
    ============================================================================================*/
    int devID = 0;
    const int sizeMult = 5;
    sMatrixSize matrix_size;

    //we start with 2^r particles. Init cude will make h_A 
    //have dimensions (#components)*numberPoints*block_size = 6*Np*2^5,
    //or an effective number of particles Np*2^5, so Np*2^5=2^r
    //implies Np=2^(r-5)
    //We're gonna have to do this in batches, with log_2(batch size)= rBATCH
    //h_A will be identical for each batch, but h_B will need
    //to be updated each batch. We can initialise cuda with 
    //identical constant parameters since the batch specifications
    //will be the same for all batches
    int np = pow(2, rBATCH-5);
    initializeCUDA(argc, argv, devID, np, matrix_size);

    
    //aight so we need this mapping to correspond to 2^r total points,
    //so with np=2^(r-5), np*2^5 should do the trick
    char type = 'c';
    Map Mcav(6 * 6, type, sBATCH);
    Mcav.Set_L(10.0f);
    Mcav.Set_Phi0(0.25f);
    Mcav.Set_Omega(1e9);
    Mcav.Set_B0(1e-3);
    Mcav.Set_v0(1.5e8);

    /*============================================================================================
    END MATRIX MAPPING INITIALISATION
    ============================================================================================*/
    /*============================================================================================
    BEGIN COMPUTATION 
    ============================================================================================*/
    Mcav.Update(beam, 0.0, type);
    float* data = Mcav.GetData();
    Mcav.ScaleMap();
    float* h_A = Mcav.GetScaled(); //this is h_A !
    assert(h_A != NULL);
    //for (unsigned int i = 0; i < 6 * Mcav.GetNp(); ++i)
    //{
    //    for (unsigned int j = 0; j < 6 * Mcav.GetNp(); ++j)
    //    {
    //        //std::cout << std::left << std::setw(8) << *(h_A + 6 * Mcav.GetNp() * i + j);
    //        std::cout << std::left << std::setw(8) << h_A[6 * Mcav.GetNp() * i + j];

    //        std::cout << "  ";
    //    }
    //    std::cout << std::endl;
    //}
    //need to get h_B which is column vector of (x1,y1,z1,vx1,vy1,vz1,x2,....,vyn,vzn)
    //inits h_B assuming the dimensions spit out by init cuda match what's been prepared by Map class
    //This is where the specification as to what batch we're doing will come in
    //only h_B is a batch dependent parameter
    //we'll assign each batches results into a second set
    //called beamOUT, so clearing h_B and h_C each batch won't be an issue;
    printf("[Matrix Multiply CUBLAS] - Starting...\n");
    for (unsigned int BATCH = 0; BATCH < nBATCH; ++BATCH) {
        if (BATCH % 5 == 0)std::cout << BATCH << "/" << nBATCH << "\r";
        const unsigned int size_B = matrix_size.uiWB * matrix_size.uiHB;
        float* h_B = new float[size_B] {0.0f};
        //assign host memeory for result per batch
        const unsigned int size_C = matrix_size.uiWC * matrix_size.uiHC;
        float* h_C = new float[size_C] {0.0f}; //this is passed by reference to matmul, so we'll be able to access the result on host device
        unsigned int cnt = 0;
        while (cnt < 6 * Mcav.GetNp()) {
            const int idx = BATCH * sBATCH + cnt % 6;
            h_B[cnt] = beam[idx].x;
            h_B[cnt + 1] = beam[idx].y;
            h_B[cnt + 2] = beam[idx].z;
            h_B[cnt + 3] = beam[idx].GB_x;
            h_B[cnt + 4] = beam[idx].GB_y;
            h_B[cnt + 5] = beam[idx].GB_z;
            cnt += 6;
        }
        //do the computation ...
        const int matrix_result = matrixMultiply(argc, argv, h_A, h_B, h_C, devID, matrix_size);
        //for (unsigned int i = 0; i < 6; ++i)
        //{
        //    for (unsigned int j = 0; j < 6 ; ++j)
        //    {
        //        //std::cout << std::left << std::setw(8) << *(h_A + 6 * Mcav.GetNp() * i + j);
        //        std::cout << std::left << std::setw(8) << h_A[6* Mcav.GetNp()* i + j];

        //        std::cout << " | ";
        //    }
        //    std::cout << "|" << h_B[i] << std::right<<std::setw(8)<<"|" << std::left <<std::setw(8) << "=? " << std::left<<std::setw(8) << h_C[i] << std::endl;
        //}
        unsigned int cnt2 = 0;
        while (cnt2 < 6 * Mcav.GetNp()) {
            const int idx = BATCH * sBATCH + cnt2 % 6;
            beamOUT[idx].x = h_C[cnt2];
            beamOUT[idx].y = h_C[cnt2 + 1];
            beamOUT[idx].z = h_C[cnt2 + 2];
            beamOUT[idx].GB_x = h_C[cnt2 + 3];
            beamOUT[idx].GB_y = h_C[cnt2 + 4];
            beamOUT[idx].GB_z = h_C[cnt2 + 5];
            cnt2 += 6;
        }
        //free host memory of A and B
        //free(h_A); //taken care of by destructor of the map
        delete[] h_B;
        delete[] h_C;
    }
    return 0;
}

//initializeCUDA_random(argc, argv, devID, sizeMult, matrix_size);
//int matrix_result = matrixMultiply_random(argc, argv, devID, matrix_size);

/*const unsigned int size_A = matrix_size.uiWA * matrix_size.uiHA;
    const unsigned int mem_size_A = sizeof(float) * size_A;
    float* h_A = (float*)malloc(mem_size_A);
    const unsigned int size_B = matrix_size.uiWB * matrix_size.uiHB;
    const unsigned int mem_size_B = sizeof(float) * size_B;
    float* h_B = (float*)malloc(mem_size_B);*/

    //srand(2006);
     //randomInit(h_A, size_A);
     //randomInit(h_B, size_B);
