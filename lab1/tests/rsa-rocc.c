//see LICENSE for license
// The following is a RISC-V program to test the functionality of the
// rsa RoCC accelerator.
// Compile with riscv-unknown-elf-gcc rsa-rocc.c
// Run with spike --extension=rsa pk a.out

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include "rsa.h"

typedef struct uint128 {
    uint64_t hi;
    uint64_t lo;
} uint128;



int main() {
/* Private-Key: (128 bit)                                                                                                                                         */
/* modulus: */
/*    00:e0:37:d3:5a:8b:16:0e:b7:f1:19:19:bf:ef:44: */
/*    09:17 */
/* publicExponent: 65537 (0x10001) */
/* privateExponent: */
/*    00:ca:b1:0c:ca:a4:43:7b:67:11:c9:77:a2:77:fe: */
/*    00:a1 */
/* prime1: 18125493163625818823 (0xfb8aafffd4b02ac7) */
/* prime2: 16442969659062640433 (0xe43129c94cf45f31) */
/* exponent1: 5189261458857000451 (0x4803f5cd8dcbfe03) */
/* exponent2: 12850891953204883393 (0xb2578a24fdb3efc1) */
/* coefficient: 10155582946292377246 (0x8cefe0e210c5a69e) */

    //DO NOT MODIFY
    uint128 modulus = {0xe037d35a8b160eb7LL,  0xf11919bfef440917LL};
    uint128 privateExp = {0x00cab10ccaa4437b67LL,  0x11c977a277fe00a1LL};
    uint64_t pubExp = 65537;
    const char plaintext[] = "Hella !";
    uint128 ciphertext;
    uint128 decrypted;
    //END DO NOT MODIFY

    //ADD TEST CASES HERE TO TEST THE CIPHERTEXT (INITIALIZE MODULUS, PLAINTEXT, PUBLIC AND PRIVATE KEY AGAIN HERE
    //AND COMPARE WITH SOME CIPHERTEXT THAT YOU CALCULATED)

    int dummy_result;

    uint64_t initCycle, duration;
    initCycle = rdcycle();

    asm volatile ("fence"); //NOTE that fences are only needed if your accelerator accesses memory
    ROCC_INSTRUCTION(0, dummy_result, modulus.lo, modulus.hi, 0);
    ROCC_INSTRUCTION(0, dummy_result, privateExp.lo, privateExp.hi, 1);
    ROCC_INSTRUCTION(0, dummy_result, &plaintext, sizeof(plaintext)/sizeof(plaintext[0]), 2);
    ROCC_INSTRUCTION(0, dummy_result, pubExp, 0, 3);
    ROCC_INSTRUCTION(0, dummy_result, &ciphertext.lo, &ciphertext.hi, 4);
    /* YOUR CODE HERE: Invoke your RSA acclerator, write the encrypted output of plaintext to ciphertext */
    asm volatile ("fence");



    //DO NOT MODIFY
    duration = rdcycle() - initCycle;
    printf("RSA Encryption took %llu cycles!\n", duration);
    initCycle = rdcycle();
    //END DO NOT MODIFY
    ROCC_INSTRUCTION(0, dummy_result, &decrypted.lo, &decrypted.hi, 5);
    /* YOUR CODE HERE: Invoke your RSA acclerator, write the decrypted output of ciphertext to decrypted */
    asm volatile ("fence");

    //DO NOT MODIFY
    duration = rdcycle() - initCycle;
    printf("RSA Decryption took %llu cycles!\n", duration);

    char *decrypted_text = (char*)(&decrypted);
    printf("decrypted=%s\n", decrypted_text);
    assert(strcmp(plaintext, decrypted_text) == 0);
}