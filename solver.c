/*
 * Copyright (c) 2018 Renato Montes
 * License: Apache License 2.0
 * https://www.apache.org/licenses/LICENSE-2.0
 *
 * Compile with:
 *     gcc -ansi -W -Wall -pedantic solver.c -o solver.exe
 *
 * Run as:
 *     ./solver.exe 3 < three.csv
 *     ./solver.exe 4 < four.csv
 *     ./solver.exe 10 < ten.csv
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Number of arguments expected in program invocation from command-line. */
static int ARG_NUMBER = 2;
/* Minimum matrix size. */
static int MIN_MATRIX_SIZE = 2;

/* Input buffer for lines in the CSV value source. */
#define BUFFSIZE 1024
/* Length of decimal values in CSV source, including negative sign. */
#define DECIMAL_LENGTH 6

/* Function prototypes. */
int     get_size(int argc, char** argv);
struct problem* get_data(int size);
void    save_constant(float* constants, int row, char line[BUFFSIZE], size_t position);
float** get_cofactors(float** original, int size);
float   get_det(float** original, int size);
void    get_minor(float** matrix, float** minor, int size, int i, int j);
float   get_det_with_cofactors(float** original, float** cofactors, int size);
void    free_matrix(float** matrix, int size);
float*  get_answers(float det, float** cofactors, float* constants, int size);
int*    get_rounded_answers(float* answers, int size);
void    print_answers(int* answers, int size);

/* Problem data entered. */
struct problem {
    float** original;
    float* constants;
};

/*
 * Check command-line arguments for a valid matrix size.
 *
 * Arguments:
 *   argc (int): number of arguments at command-line invocation
 *   argv (char**): string arguments at command-line invocation
 *
 * Return:
 *   int: size entered as a command-line argument
 *
 * Post-condition:
 *   the returned matrix size is valid
 */
int get_size(int argc, char** argv) {
    int size;
    if (argc != ARG_NUMBER) {
        fprintf(stderr, "Usage: speedy [size] < [datafile]\n");
        exit(0);
    }
    if ((size = atoi(argv[1])) < MIN_MATRIX_SIZE) {
        fprintf(stderr, "Error: matrix size must be an integer greater"\
                        " than or equal to 2\n");
        exit(0);
    };
    return size;
}

/* 
 * Obtain CSV data from stdin.
 *
 * Arguments:
 *   size (const int): size of the square matrix
 *
 * Return:
 *   struct problem: data obtained from CSV input
 */
struct problem* get_data(int size) {
    int row, column;
    char line[BUFFSIZE];
    size_t position;
    size_t counter;
    char linechar;
    float linefloat;
    char substring[DECIMAL_LENGTH];
    struct problem* data = malloc(sizeof(struct problem));
    float** matrix = malloc(sizeof(float*) * size);
    float* constants = malloc(sizeof(float) * size);

    for (row = 0; row < size; row++) {
        if (!fgets(line, BUFFSIZE, stdin)) {
            fprintf(stderr, "Error: could not read a line\n");
            exit(0);
        }
        matrix[row] = malloc(sizeof(float) * size);
        position = 0;
        for (column = 0; column < size; column++) {
            counter = 0;
            while ((linechar = line[position]) != ',') {
                substring[counter] = linechar;
                counter++;
                position++;
            }
            substring[counter] = '\0';
            sscanf(substring, "%f", &linefloat);
            matrix[row][column] = linefloat;
            position++;
        }
        counter = 0;
        while ((linechar = line[position]) != '\n' && linechar != '\0') {
            substring[counter] = linechar;
            counter++;
            position++;
        }
        substring[counter] = '\0';
        sscanf(substring, "%f", &linefloat);
        constants[row] = linefloat;
    }
    data->original = matrix;
    data->constants = constants;
    return data;
}

/*
 * Calculate cofactor matrix from original matrix.
 *
 * Arguments:
 *   original (float**): original matrix
 *   size (int): size of the square matrix
 *
 * Return:
 *   float**: cofactor matrix, including correct positive/negative signs
 */
float** get_cofactors(float** original, int size) {
    int row, col;
    int sign;
    float** cofactors = malloc(sizeof(float*) * size);
    float** topminor = malloc(sizeof(float*) * (size - 1));

    for (row = 0; row < size; row++) {
        cofactors[row] = malloc(sizeof(float) * size);
    }

    for (row = 0; row < (size - 1); row++) {
        topminor[row] = malloc(sizeof(float) * (size - 1));
    }

    for (row = 0; row < size; row++) {
        sign = (row % 2 == 0)? 1 : -1;
        for (col = 0; col < size; col++) {
            get_minor(original, topminor, size, row, col);
            cofactors[row][col] = sign * get_det(topminor, size - 1);
            sign = -sign;
        }
    }

    for (row = 0; row < (size - 1); row++) {
        free(topminor[row]);
    }
    free(topminor);

    return cofactors;
}

/*
 * Calculate determinant of matrix in a recursive manner.
 *
 * Arguments:
 *   matrix (float**): the square matrix to obtain the determinant of
 *   size (int): size of the square matrix
 *
 * Return:
 *   float: the determinant of the input matrix
 */
float get_det(float** matrix, int size) {
    float det;
    int j;
    int sign = 1;
    float** minor = malloc(sizeof(float*) * (size - 1));
    for (j = 0; j < size - 1; j++) {
        minor[j] = malloc(sizeof(float) * (size - 1));
    }

    if (size == 2) {
        return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }

    for (j = 0; j < size; j++)
    {
        get_minor(matrix, minor, size, 0, j);
        det += sign * matrix[0][j] * get_det(minor, size - 1); 
        sign = -sign;
    }

    for (j = 0; j < size - 1; j++) {
        free(minor[j]);
    }
    free(minor);

    return det;
}

/*
 * Obtain minor based on specific element and modify in-place.
 *
 * Arguments:
 *   matrix (float**): the matrix to examine
 *   minor (float**): the matrix to store the minor in
 *   size (int): the size of the matrix to examine
 *   i (int): the row of the element of reference to discriminate
 *   j (int): the column of the element of reference to discriminate
 */
void get_minor(float** matrix,
               float** minor,
               int size,
               int i,
               int j) {
    int matrixrow;
    int matrixcol;
    int minorrow = 0;
    int minorcol = 0;

    for (matrixrow = 0; matrixrow < size; matrixrow++) {
        for (matrixcol = 0; matrixcol < size; matrixcol++) {
            if (matrixrow != i && matrixcol != j) {
                minor[minorrow][minorcol] = matrix[matrixrow][matrixcol];
                minorcol++;
                if (minorcol == size - 1) {
                    minorcol = 0;
                    minorrow++;
                }
            }
        }
    }
}

/*
 * Get a determinant quickly by making use of a cofactor matrix.
 * 
 * Arguments:
 *   original (float**): the original matrix
 *   cofactors (float**): the cofactor matrix
 *   size (int): the size of the square matrices
 *
 * Return:
 *   float: the calculated determinant
 */
float get_det_with_cofactors(float** original,
                             float** cofactors,
                             int     size) {
    int col;
    float det = 0.0;
    for (col = 0; col < size; col++) {
        det += original[0][col] * cofactors[0][col];
    }
    return det;
}

/*
 * Calculate answers to the AX = B system of linear equations
 *
 * Arguments:
 *   det (float): the determinant of matrix A
 *   cofactors (float**): the cofactors of matrix A
 *   constants (float*): the right-hand constants (B)
 *   size (int): the size of the square matrix A
 *
 * Return:
 *   float*: answers to the linear equations, as floats
 */
float* get_answers(float   det,
                   float** cofactors,
                   float*  constants,
                   int     size) {
    int i, j;
    float* answers = malloc(sizeof(float) * size);
    for (i = 0; i < size; i++) {
        answers[i] = 0;
        for (j = 0; j < size; j++) {
            /* Inverse matrix built on the fly. */
            answers[i] += constants[j] * cofactors[j][i];
        }
        answers[i] = answers[i] / det;
    }
    return answers;
}

/*
 * Round float answers into ints.
 *
 * Arguments:
 *   answers (float*): answers obtained as floats
 *   size (int): number of answers obtained
 *
 * Return:
 *   int*: answers obtained, now rounded to be ints
 *
 * Pre-condition:
 *   All answers for the input data are intended to be ints. This function
 *   approximates all received floats to their closest int value.
 */
int* get_rounded_answers(float* answers, int size) {
    int i;
    float cur;
    int* final = malloc(sizeof(int) * size);
    for (i = 0; i < size; i++) {
        cur = answers[i];
        final[i] = (cur < 0.0)? (int)(cur - 0.5) : (int)(cur + 0.5) ;
    }
    return final;
}

/*
 * Print int answers to the command-line.
 *
 * Arguments:
 *   answers (int*): answers to the problem data, as ints
 *   size (int): number of answers to the problem data (number of equations)
 */
void print_answers(int* answers, int size) {
    int i;
    printf("%d", answers[0]);
    for (i = 1; i < size; i++) {
        printf(",%d", answers[i]);
    }
}

/*
 * Free the allocated memory of a matrix.
 *
 * Arguments:
 *   matrix (float**): matrix to be freed
 *   size (int): size of the matrix to be freed
 */
void free_matrix(float** matrix, int size) {
    int row;
    for (row = 0; row < size; row++) {
        free(matrix[row]);
    }
    free(matrix);
}

/*
 * Program entry point.
 *
 * Arguments:
 *   argc (int): number of command-line arguments at program invocation
 *   argv (char**): strings of command-line arguments at program invocation
 *
 * Return:
 *   int: program runtime exit code
 */
int main(int argc, char** argv) {
    int size;
    struct problem data;
    float** cofactors;
    float det;
    float* answers;

    /* Read matrix size from command-line argument. */
    size = get_size(argc, argv);
    /* Read data from input. */
    data = *(get_data(size));
    /* Obtain cofactor matrix. */
    cofactors = get_cofactors(data.original, size);
    /* Obtain determinant using the cofactor matrix. */
    det = get_det_with_cofactors(data.original, cofactors, size);
    /* Calculate answers as floats. */
    answers = get_answers(det, cofactors, data.constants, size);
    /* Print answers after rounding them to the closest int values. */
    print_answers(get_rounded_answers(answers, size), size);
    /* Free last matrices allocated in memory. */
    free_matrix(data.original, size);
    free_matrix(cofactors, size);

    return 0;
}


