//============================================================
//
//      Type:        FWH tools file
//
//      Author:      Tommaso Bellosta on 14/04/26.
//                   Dipartimento di Scienze e Tecnologie Aerospaziali
//                   Politecnico di Milano
//                   Via La Masa 34, 20156 Milano, ITALY
//                   e-mail: tommaso.bellosta@polimi.it
//
//      Copyright:   2026, authors above and the FWH contributors.
//                   This software is distributed under the MIT license, see LICENSE.txt
//
//============================================================
#include "geometry.h"

/** a function. Parameters are passed by reference
 * i.e. the function operates on the actual elements it gets passed,
 * not on copies. If you modofy them inside the function, they will
 * be changes also on the caller side.**/
double multiply(double &a, double &b) {
    double c = 3.141592653589793;


    a = a * 10;
    c = a * b;
    return c;
}

/** a function. Parameters are passed "by value"
 * i.e. the function operates on copies of the inputs it gets passed.
 * If you modify them inside the function, they will
 * NOT change on the caller side.**/
double multiply2(double a, double b) {
    double c = 3.141592653589793;


    a = a * 10;
    c = a * b;
    return c;
}

/** a function. Parameters are passed by const reference
 * i.e. the function operates on the actual elements it gets passed,
 * not on copies. Since inputs are const references, you will not be
 * able to modify them in the function body.**/
double multiply3(const double &a, const double &b) {
    double c = 3.141592653589793;

    // this does not compile!
    // a = a*10;
    c = a * b;
    return c;
}

int main() {
    /** define a variable of type double.
     * You can also initialize a variable at creation:
     * -- double d = 0;
     * -- double d{0};
     * You can also define multiple variables in the same line:
     * -- double d,f=0,g,h=2.2;
     **/
    double d;

    /** standard operation on a double variable **/
    d = 4.1;
    d = d * 3;
    d += 3;

    /** You can allocate a contiguous block of memory of size
     * known at compile time (i.e. a constant) on the stack. This
     * is called an array. **/
    double dArray[5];

    /** fill the elements in the array with 1..5
     * the square brackets give you access to an
     * element of the array. Make sure not to access elements
     * outside the memory allocated by the array. By doing so you will
     * corrupt the memory and the program will give you UD
     * (undefined behaviour, meaning there is no way to predict what the
     * program will do. Hopefully it will crash..)
     **/
    for (int i = 0; i < 5; i++) {
        dArray[i] = i + 1;
    }

    /** The following line creates a pointer of type double.
     * A pointer holds a memory address. You can point it
     * to an already allocated region of memory, and use it to handle
     * the data store at that address. **/
    double *dptr; // this is an empty pointer. it does not point to anything

    /** assign to the pointer the address in memory where variable
     * d is stored. You can get the memory address of a variable by
     * using "&" in front of the variable name. **/
    dptr = &d;

    /** Now, since ptrs stores the address of variable d, the following
     * two line are equivalent. When a pointer stores the address in memory of
     * a variable, you can get the variable it points to by using
     * a "*" in front of the pointer (i.e. *dptr gives you d) **/
    d = 3;
    *dptr = 3;


    /** We can use pointers also to point to arrays. In fact, arrays are just
     * pointers to the beginning of the memory allocation for that array,
     * i.e. dArray is a pointer of type double storing the address in memory
     * of the contiguous block of memory of 5 doubles. **/

    /** the following two lines are equivalent.
     *  We are pointing dptr to the beginning of the memory addressed by
     *  dArray. We can do that by taking the address of the first element
     *  in dArray (&dArray[0]), or by taking the pointer to the array
     *  allocation (i.e. dArray is a pointer to a stack allocated block of memory)**/
    dptr = &dArray[0];
    dptr = dArray;

    /** again, equivalent **/
    dArray[4] = 76.25;
    dptr[4] = 76.25;

    /** The size of arrays living on the stack must be known at
     * compile time, i.e. you must know the size in advance. We want
     * to fix this, i.e. we want to be able to define arrays at runtime,
     * meaning with a size that depends on the code execution (or inputs).
     * We can do so by declaring our arrays on the heap (another region
     * of memory). We first create a pointer, then we assign to it a
     * contiguous region of memory on the heap. The syntax is:
     * -- type* = new type[size]
     **/
    double *f;
    f = new double[84]; //creates a heap array of 84 doubles.
    /** can also do it in a single line:
     * -- double *f = new double[84];
     * -- auto   *f = new double[84];
     **/

    /** You can use it as you do stack arrays **/
    f[64] = 1.1;
    for (int i = 0; i < 84; i++) {
        f[i] = 0;
    }
    /** A while loop **/
    int i = 0;
    while (i < 84) {
        f[i] = 0;
        i++;
    }

    /** heap allocated memory need to be deallocated manually  **/
    delete [] f;

    /** References are aliases of already existing variables.**/
    double &ref = d;
    /** the following lines are equivalent **/
    d = 5;
    ref = 5;

    /** a const reference aliases an existing varialbe but cannot
     * overwrite its content **/
    const double &cref = d;
    double g = 6 * cref;
    // cref = 3.14; // does not compile


    /** You can also create single variables on the
     * stack, not just blocks of variables: **/
    double *single = new double;
    delete single;

    /** call a function **/
    double a, b;
    double c = multiply(a, b);


    /** Objects are instances of a class. You can define and
     * use objects as you did with simple types (int, double..).
     * In order to use our Vec3 class we need to tell the compiler where
     * we declared it. The compiler only needs the function signature,
     * i.e. the return type and the types of the inputs (this is why we
     * are able to separate functions declaration in .h files and definitions
     * in .cpp files). To give the compiler this info we need to include the
     * declaration of the class with the following line (usually placed
     * at the beginning of the file):
     * -- #include "geometry.h"
     **/

    /** create a Vect3 object. This calls the first constructor
     * that sets all 3 components to zero **/
    Vect3 pos;
    /** create another Vect3 object. This calls the secon constructor
     * that sets custom values to the components **/
    Vect3 pos2(1, 1, 1);

    /** We defined the [] symbols for our Vect3 class to return
     * the components of the 3D vector. It either returns a reference
     * to the i-th component or a const reference, depending on the
     * type of the operation we are doing.
     **/
    pos[1] = 0; //the [] operator returns a reference
    double& y = pos[2]; //the [] operator returns a reference

    double val = pos[0]; // the [] operator returns a const reference
    double x = 3.14 * pos[3]; //  the [] operator returns a const reference

    /** We also defined the + operation between the Vect3 objects **/
    Vect3 pos3 = pos + pos2;
}
