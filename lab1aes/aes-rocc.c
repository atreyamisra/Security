//see LICENSE for license
// The following is a RISC-V program to test the functionality of the
// aes RoCC accelerator.
// Compile with riscv-unknown-elf-gcc aes-rocc.c
// Run with spike --extension=aes pk a.out

#include <assert.h>
#include <stdio.h>
#include <stdint.h>
#include "aes.h"
struct AES_ctx
{
    uint8_t RoundKey[240];
    uint8_t Iv[16];
};
typedef uint8_t state_t[4][4];
int main() {
    //DO NOT MODIFY
    unsigned char enc_buf[128];
    unsigned char plaintext[1][32] = {
	    {0x6b,0xc1,0xbe,0xe2,0x2e,0x40,0x9f,0x96,0xe9,0x3d,0x7e,0x11,0x73,0x93,0x17,0x2a,0xae,0x2d,0x8a,0x57,0x1e,0x03,0xac,0x9c,0x9e,0xb7,0x6f,0xac,0x45,0xaf,0x8e,0x51}
    };
    unsigned char ciphertext[1][32] = {
	    {0x60,0x1e,0xc3,0x13,0x77,0x57,0x89,0xa5,0xb7,0xa7,0xf5,0x04,0xbb,0xf3,0xd2,0x28,0xf4,0x43,0xe3,0xca,0x4d,0x62,0xb5,0x9a,0xca,0x84,0xe9,0x90,0xca,0xca,0xf5,0xc5}
    };
    unsigned char iv[1][16] = {
	    {0xf0,0xf1,0xf2,0xf3,0xf4,0xf5,0xf6,0xf7,0xf8,0xf9,0xfa,0xfb,0xfc,0xfd,0xfe,0xff},
    };
    unsigned char key[1][32] = {
	    {0x60,0x3d,0xeb,0x10,0x15,0xca,0x71,0xbe,0x2b,0x73,0xae,0xf0,0x85,0x7d,0x77,0x81,0x1f,0x35,0x2c,0x07,0x3b,0x61,0x08,0xd7,0x2d,0x98,0x10,0xa3,0x09,0x14,0xdf,0xf4}
    };

    unsigned char decrypted_text[32];
    //END DO NOT MODIFY

    int dummy_result;


    unsigned long long initCycle, duration;
    initCycle = rdcycle();
    asm volatile ("fence"); //NOTE: fences are only needed if the accelerator is accessing memory
    ROCC_INSTRUCTION(0, dummy_result, &enc_buf, &plaintext[0], 0);
    ROCC_INSTRUCTION(0, dummy_result, &ciphertext[0], &iv[0], 1);
    ROCC_INSTRUCTION(0, dummy_result, &key[0], &decrypted_text, 2);
    ROCC_INSTRUCTION(0, dummy_result, 0, 0, 3);
    //printf("%d", n);

    //YOUR CODE HERE: Invoke your AES acclerator, write the encrypted output of plaintext to enc_buf
    asm volatile ("fence");

    //DO NOT MODIFY
    duration = rdcycle() - initCycle;
    printf("AES Encryption took %llu cycles!\n", duration);
    initCycle = rdcycle();
    //END DO NOT MODIFY
    asm volatile ("fence");
    ROCC_INSTRUCTION(0, dummy_result, 0, 0, 4);
    //YOUR CODE HERE: Invoke your AES acclerator, write the decrypted output of enc_buf to decrypted_text
    asm volatile ("fence");

    //DO NOT MODIFY
    duration = rdcycle() - initCycle;
    printf("AES Decryption took %llu cycles!\n", duration);

    // Check result
    assert(memcmp(enc_buf, ciphertext[0], 32) == 0);
    assert(memcmp(decrypted_text, plaintext[0], 32) == 0);
    //END DO NOT MODIFY
    return 0;
}