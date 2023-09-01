#include <iostream>
#include <cmath>
using namespace std;

int binary2decimalfun(int binary[], int size) {
    int decimal = 0;
    for (int i = size - 1; i >= 0; i--) {
        decimal += binary[i] * pow(2, size - i - 1);
    }
    return decimal;
}

void binary2decimalobj(int binary[], int size, int& decimal) {
    for (int i = size - 1; i >= 0; i--) {
        decimal += binary[i] * pow(2, size - i - 1);
    }
}

void decimal2binary(int decimal, int objekt[]) {
    const int length = ceil(log2(decimal+1));
    for (int i; i >= 0; i--) {
        if (decimal >= pow(2, i))  {
            objekt[length - i - 1] = 1;
            decimal -= pow(2,i);
        }
    }    
}

int main() {
    int size = 8;
    int decimal = 0;
    
    int binary[] = {1, 0, 1, 0, 1, 0, 1, 0};

    binary2decimalobj(binary, size, decimal);

    cout << decimal << endl;
    int objekt[] = {0,0,0,0,0,0,0,0};
    decimal2binary(decimal, objekt);

    for(int i=0; i < 8; i++) {

        cout << objekt[i] << " ";
    }
    return 0;
}