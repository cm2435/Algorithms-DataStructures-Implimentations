#include <fstream>
#include <iostream>
#include <math.h>
#include <stdlib.h>

//define pi for the h1,h2 functions
long double pi = 3.14159265358979323846264338327;
long double pi_2 = pi / 2;

//Using namespace std for brevvity since only small file, would be bad convention for larger projects
using namespace std;

//define imaginary number functionality
typedef struct complex_number
{
    long double real;
    long double imaginary;
    long double time;
};

//define struct for reading in h3
typedef struct h3
{
    size_t N;
    long double time;
    complex_number complex_reading;
} h3;

//End of struct definitions
//##############################################################################

//define "pass" function analagous to python for use later on
void pass() {}

//basic complex number logic functions to compute the mod of complex numbers
double complex_mod(complex_number z)
{
    return sqrt(z.real * z.real + z.imaginary * z.imaginary);
}

double mod_array(complex_number *z)
{
    return complex_mod(*z);
}

//function for comparing the amplitude of mod of complex_numbers for using Qsort
int CompareComplex(void *a, void *b)
{
    double mod_a = mod_array((complex_number *)a);
    double mod_b = mod_array((complex_number *)b);

    int x = 0;

    if (mod_a > mod_b)
    {
        x = 1;
    }
    else if (mod_b > mod_a)
    {
        x = -1;
    }
    else
    {
        x = 0;
    }

    return x;
}

//extention of compare complex for quicksort with h3
int Compare_h3(const void *a, const void *b)
{
    //covert a,b into h3 structs
    const h3 *a_2 = (h3 *)a;
    const h3 *b_2 = (h3 *)b;

    return CompareComplex((void *)&a_2->complex_reading, (void *)&b_2->complex_reading);
}

//check position in array helper for later on
int checkposition(size_t *ls, size_t sz, size_t x)
{
    size_t i;
    for (i = 0; i < sz; ++i)
    {
        if (ls[i] == x)
        {
            return 1;
        }
    }
    return 0;
}

//define h1 function
complex_number h1_Complex(const double t)
{
    complex_number b;
    b.real = cos(t) + cos(5 * t);
    b.imaginary = sin(t) + sin(5 * t);
    b.time = t;

    //cout << "output of the function H1 is :  " << b.real << b.imaginary << "i" << "\n\n";
    return b;
}

//define h2 function
complex_number h2_Complex(const double t)
{
    complex_number c;

    c.real = exp(pow(t - pi, 2) / 2);
    c.imaginary = 0;
    c.time = t;

    return c;
}

//define a sampling function analagous to pythons numpy.linspace
double *linspace(double start, double end, size_t N)
{
    // Allocate memory for result, calculate increment size for N values in given range
    double *arr = (double *)malloc(N * sizeof(double));

    // Calculate increment for N evenly spaced values from start -> end
    // 'val' keeps track of current value
    double increment = (double)(end - start) / N;
    double val = start;

    // Add values to arr, incrementing val every time
    size_t i;
    for (i = 0; i < N; i++)
    {
        arr[i] = val;
        val += increment;
    }
    return arr;
}

//function to compute the discrete fourier transform
complex_number *DFT(const complex_number *sample, const size_t N)
{
    /*
    params- sample: array of structs, was initially a vector but rewrote to not use c++ concepts
            N: the number of elements in the array

    returns- array of DFT transformed structs of len = input

    */

    complex_number *arr = (complex_number *)malloc(N * sizeof(complex_number));
    complex_number H_k, H_j;
    double exponent, exponent_j;
    size_t i, j;

    for (i = 0; i < N; i++)
    {
        H_k.imaginary = 0.;
        H_k.real = 0.;
        //calc to avoid doing in nested loop
        exponent = (-2. * pi * i / N);

        for (j = 0; j < N; ++j)
        {
            exponent_j = exponent * j;

            //Get j'th element from sample
            H_j = sample[j];

            H_k.real += (H_j.real * cos(exponent_j) + H_j.imaginary * sin(exponent_j));
            H_k.imaginary += (H_j.imaginary * cos(exponent_j) - H_j.real * sin(exponent_j));
        }

        arr[i] = H_k;
    }
    return arr;
}

//perform the inverse fourier transform
complex_number *IFT(complex_number *sample, size_t N, size_t *skip_n, size_t size)
{
    /*
    params- sample: array of pointers to structs, was initially a vector but rewrote to not use c++ concepts
            N: the number of elements in the array
            skip_N: the number of initial samples from the array to skip.
            size: size of index to be used for position checking

    returns- array of IFT transformed structs of len = input

    */

    complex_number H_k, H_j;
    double exponent, exponent_j;
    complex_number *arr = (complex_number *)malloc(N * sizeof(complex_number));
    size_t i, j;

    for (j = 0; j < N; j++)
    {
        H_k.imaginary = 0.;
        H_k.real = 0.;

        //precalc again outside nested loop
        exponent = (2. * pi * j / N);

        for (i = 0; i < N; i++)
        {
            //Check the indexed position to allow for skipping downstream
            if (checkposition(skip_n, size, i))
            {
                pass();
            } else {
                exponent_j = exponent * i;

                //retrieve jth sample
                H_j = sample[i];

                H_k.real += ((H_j.real * cos(exponent_j) + H_j.imaginary * sin(exponent_j)));
                H_k.imaginary += ((H_j.imaginary * cos(exponent_j) - H_j.real * sin(exponent_j)));
            }
        }
        H_k.real = H_k.real / N;
        H_k.imaginary = H_k.imaginary / N;

        arr[j] = H_k;
    }

    return arr;
}

//End of general functions, now question implimentations
//####################################################################################################################################################

//Begin direct implimentations for the various parts

complex_number **Q_3b(double a, const double b, size_t N)
{
    double *sample = linspace(a, b, N);
    complex_number **results = (complex_number **)malloc(2 * sizeof(complex_number *));

    complex_number *output_h1 = (complex_number *)malloc(N * sizeof(complex_number));
    complex_number *output_h2 = (complex_number *)malloc(N * sizeof(complex_number));

    //copy output of the h1,h2 functions sampled from 0-2pi to new array
    for (int i = 0; i < N; ++i)
    {
        output_h1[i] = h1_Complex(sample[i]);
        output_h2[i] = h2_Complex(sample[i]);
    }

    //Output to files
    ofstream file1, file2;
    file1.open("complex_sample_h1.txt");
    file1 << "time, real part, imaginary part" << endl;

    file2.open("complex_sample_h2.txt");
    file2 << "time, real part, imaginary part" << endl;

    for (int i = 0; i < 100; ++i)
    {
        file1 << output_h1[i].time << ",  (" << output_h1[i].real << "+" << output_h1[i].imaginary << "i)" << endl;
        file2 << output_h2[i].time << ",  (" << output_h2[i].real << "+" << output_h2[i].imaginary << "i)" << endl;
    }

    file1.close();
    file2.close();

    //save to results array to feed into next question
    results[0] = output_h1;
    results[1] = output_h2;

    return results;
}

complex_number **Q_3d_E(double *times, complex_number **samples, size_t N)
{

    size_t i;
    complex_number **results = (complex_number **)malloc(2 * (sizeof(complex_number *)));

    //Take the DFT of the input samples
    complex_number *H1 = DFT(samples[0], N);
    complex_number *H2 = DFT(samples[1], N);

    //Assign to dummy array
    results[0] = H1;
    results[1] = H2;

    //printing out data
    cout << "printing H1:" << endl;
    cout << "time, H1 real part, H1 imaginary" << endl;

    for (i = 0; i < N; ++i)
    {
        cout << times[i] << ", (" << H1[i].real << "+" << H1[i].imaginary << "i)" << endl;
    }

    cout << "printing H2:" << endl;
    cout << "H2 time, H2 real part, H2 imaginary part" << endl;

    for (i = 0; i < N; ++i)
    {
        cout << times[i] << ", (" << H2[i].real << "+" << H2[i].imaginary << "i)" << endl;
    }

    //return the values of H1, H2 for next task
    return results;
}

complex_number **Q_3F(complex_number **FT_variables, size_t N)
{
    // Unpack H1 & H2 from samples.
    complex_number *H1 = FT_variables[0];
    complex_number *H2 = FT_variables[1];

    complex_number **arr = (complex_number **)malloc(2 * sizeof(complex_number *));

    //define the skiprows--
    size_t* H1_skip_index = (size_t*)malloc(sizeof(size_t));
    H1_skip_index[0] = 1;
    size_t* H2_skip_index = (size_t*)malloc(sizeof(size_t));
    H2_skip_index[0] = 0;

    //Compute the IFT for the given inputs
    complex_number *h1_IFT = IFT(H1, N, H1_skip_index, 1);
    complex_number *h2_IFT = IFT(H2, N, H2_skip_index, 1);

    // Pack h1_prime & h2_prime into results so it can be returned for q_3f
    arr[0] = h1_IFT;
    arr[1] = h2_IFT;

    free(H1_skip_index);
    free(H2_skip_index);

    return arr;
}

void Q_3g(complex_number **data, double *time)
{
    //write out the output of 3f to text files for downstream plotting
    ofstream file3, file4;
    file3.open("inverse_1.txt");
    file3 << "times,real,imag" << endl;

    file4.open("inverse_2.txt");
    file4 << "times,real,imag" << endl;

    for (int i = 0; i < 100; ++i)
    {
        file3 << time[i] << "," << data[0][i].real << "," << data[0][i].imaginary << endl;
        file4 << time[i] << "," << data[1][i].real << "," << data[1][i].imaginary << endl;
    }

    file3.close();
    file4.close();
}

h3 *Q_3H(const char *file, size_t N)
{
    h3 *measurements = (h3 *)malloc(N * sizeof(h3));
    FILE *fp = fopen(file, "r");
    int index;
    complex_number z;
    size_t size = 0;
    double real_part = 0., imaginary_part = 0., time;

    while (fscanf(fp, "%d, %lf, %lf, %lf", &index, &time, &real_part, &imaginary_part) != EOF && size < N)
    {
        z.real = real_part;
        z.imaginary = imaginary_part;

        measurements[size].N = index;
        measurements[size].time = time;
        measurements[size].complex_reading = z;

        size++;
    }

    fclose(fp);
    return measurements;
}

h3 *Q_3I(h3 *data, const size_t N)
{
    complex_number *complex_arr = (complex_number *)malloc(N * sizeof(h3));
    h3 *results = (h3 *)malloc(N * sizeof(h3));
    size_t i;

    for (i = 0; i < N; ++i)
    {
        complex_arr[i] = data[i].complex_reading;
    }

    complex_number *DFT_sample = DFT(complex_arr, N);

    //copy transformed data into new array
    for (i; i < N; ++i)
    {
        results[i].N = data[i].N;
        results[i].time = data[i].time;
        results[i].complex_reading = DFT_sample[i];
    }

    //Free up not data not being returned
    free(complex_arr);
    free(DFT_sample);

    return results;
}

void Q_3K(h3 *samples, size_t N)
{
    complex_number *result = (complex_number *)malloc(N * sizeof(complex_number));

    h3 *samples_sorted = (h3 *)malloc(N * sizeof(h3));

    //Assign all but the first 4 values of the array to a list to skip
    size_t i, skip_n[196];

    //copy the sample data to an array to be sorted
    for (i = 0; i < N; ++i)
    {
        samples_sorted[i] = samples[i];
    }

    //perform the quicksort on the copied data
    qsort(samples_sorted, N, sizeof(h3), Compare_h3);

    //copy sorted data index to a skip index
    for (i = 0; i < N - 4; ++i)
    {
        skip_n[i] = samples_sorted[i].N;
    }

    //create array of complex values from original sample
    complex_number *complex_arr = (complex_number *)malloc(N * sizeof(complex_number));

    //Copy to a array to use to compute the IFT
    for (i = 0; i < N; ++i)
    {
        complex_arr[i] = samples[i].complex_reading;
    }

    //pass through the IFT function to get new array
    result = IFT(complex_arr, N, skip_n, 196);
    //Write to file 5
    ofstream file5;
    file5.open("inverse_3.txt");

    file5 << "N, time, real part, imaginary part" << endl;
    for (int i = 0; i < N; ++i)
    {
        file5 << skip_n[i] << "," << result[i].time << "," << result[i].real << "," << result[i].imaginary << endl;
    }

    file5.close();
}

//main function routine
int main()
{
    //define the sample sizes for the sampling
    size_t N = 100;

    //define an array of times to be plugging into FT and IFT
    double *times = linspace(0., 2 * pi, N);

    //Get an array of complex numbers for an array of times
    complex_number **h1_h2 = Q_3b(0, 2 * pi, N);

    //Define the values from the DFT
    complex_number **DFT_h1_h2 = Q_3d_E(times, h1_h2, N);

    //Take the IFT of h1,h2
    complex_number **IFT_h1_h2 = Q_3F(DFT_h1_h2, N);

    //Write out the IFT data to two seperate text files
    Q_3g(IFT_h1_h2, times);

    //Define next sample param
    N = 200;

    //Load in data from the h3 file for processing
    h3 *measure_data = Q_3H("data.txt", N);

    //Take the DFT of the h3 sample
    h3 *measure_data_FT = Q_3I(measure_data, N);

    //Take the IFT of the DFT data from above
    Q_3K(measure_data_FT, N);

    //Free the excess data
    free(times);
    free(DFT_h1_h2[0]);
    free(DFT_h1_h2[1]);
    free(IFT_h1_h2);

    return 0;
}
