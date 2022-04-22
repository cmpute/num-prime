# Performance comparison

The following benchmark results are produced by criterion (the code is under [`bench/`](./bench/)) on my Laptop (Ryzen 4700U).

```log
primality check (u64)/num-prime (this crate)                                                                            
                        time:   [750.96 us 756.01 us 763.51 us]
primality check (u64)/number-theory                                                                             
                        time:   [1.1550 ms 1.1586 ms 1.1624 ms]
primality check (u64)/num-primes                                                                            
                        time:   [429.79 ms 431.64 ms 434.03 ms]
primality check (u64)/glass_pumpkin                                                                            
                        time:   [108.30 ms 108.85 ms 109.42 ms]
primality check (u64)/primal-check                                                                            
                        time:   [17.693 ms 17.932 ms 18.215 ms]

primality check (u256)/num-prime (this crate)                                                                            
                        time:   [828.72 us 848.85 us 871.09 us]
primality check (u256)/num-primes                                                                             
                        time:   [1.0797 ms 1.0869 ms 1.0960 ms]
primality check (u256)/glass_pumpkin                                                                            
                        time:   [332.27 us 337.31 us 344.24 us]

safe primality check (u256)/num-prime (this crate)                                                                            
                        time:   [462.10 us 463.54 us 465.05 us]
safe primality check (u256)/num-primes                                                                            
                        time:   [796.00 us 798.52 us 801.12 us]
safe primality check (u256)/glass_pumpkin                                                                            
                        time:   [319.49 us 320.79 us 322.16 us]
safe primality check (u256)/glass_pumpkin (BPSW)                                                                            
                        time:   [315.92 us 317.33 us 318.77 us]

primality check (u2048)/num-prime (this crate)                                                                            
                        time:   [26.156 ms 26.424 ms 26.757 ms]
primality check (u2048)/num-primes                                                                            
                        time:   [22.424 ms 22.690 ms 23.073 ms]
primality check (u2048)/glass_pumpkin                                                                            
                        time:   [5.7849 ms 5.8078 ms 5.8308 ms]
primality check (u2048)/glass_pumpkin (BPSW)                                                                            
                        time:   [5.8614 ms 5.9120 ms 5.9677 ms]

safe primality check (u2048)/num-prime (this crate)                                                                            
                        time:   [15.183 ms 15.269 ms 15.358 ms]
safe primality check (u2048)/num-primes                                                                            
                        time:   [46.018 ms 47.346 ms 48.770 ms]
safe primality check (u2048)/glass_pumpkin                                                                            
                        time:   [5.7260 ms 5.7565 ms 5.7890 ms]
safe primality check (u2048)/glass_pumpkin (BPSW)                                                                            
                        time:   [5.6429 ms 5.6716 ms 5.7022 ms]

factorize (u64)/num-prime (this crate)                                                                            
                        time:   [9.8809 ms 9.9713 ms 10.071 ms]
factorize (u64)/number-theory                                                                            
                        time:   [38.202 ms 38.307 ms 38.417 ms]

prime generation (256 bits)/num-prime (this crate)                                                                            
                        time:   [1.5424 ms 1.5778 ms 1.6093 ms]
prime generation (256 bits)/num-primes                                                                           
                        time:   [6.0673 ms 6.4005 ms 6.7287 ms]
prime generation (256 bits)/glass_pumpkin                                                                            
                        time:   [2.3650 ms 2.4890 ms 2.6202 ms]

safe prime generation (256 bits)/num-prime (this crate)                        
                        time:   [116.29 ms 168.60 ms 232.87 ms]
safe prime generation (256 bits)/num-primes                                                                          
                        time:   [1.1444 s 1.6604 s 2.1542 s]
safe prime generation (256 bits)/glass_pumpkin                                                                          
                        time:   [213.55 ms 432.20 ms 677.55 ms]
```