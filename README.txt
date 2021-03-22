changes made
file name: MIMO_OFDM_TX.m

i) snr is provided by Rx code where Tx code is called every monte carlo iteration and mod order iteration.
ii) variable "mul_tx"(used in line 94 and 140) is introduced to generate 2 stream of TX(TX A and Tx B as in MIMO case) or a single stream (both the TX A and B are same as in SIMO case i.e., Q7.2 and Q7.3)
iii) rx_vec_air (RX receive stream generation) in line 141 to 145 are changed in sign of addition to simulate the out of phase addition for Q7.3 case

file name: MIMO_OFDM_RX_test_SIMO_case.m

i)iter in line 4 decides the monte carlo iteration number for each snr and mod orders
ii)snr_iter in line 5 decides the total iteration of snr steps from starting value(starting value is given in line 6)
iii)mul_tx is declared to be zero for simo case
iv)h1, h2, h3 and h4(in line 126,127,128,129) are channel matrix for 4x1 case(Q7.2)
v)h1_2 and h2_2(in line 137 and 138) are channel matrix for 2x2 systems simo as in Q7.3
vi) MRC is used
vii)rx_data_final_1 and rx_data_final_2 are decoded rx data for 4x1(Q7.2) and 2x2(Q7.3) case
viii) MIMO_OFDM_TX is called each time for data generation in monte carlo simiulation
ix)ber is BER for Q7.2
x)ber1 is BER for Q7.3

file name: MIMO_OFDM_RX_test_2_TX_case.m

i) All the parameters are similar to defines above with mul_tx is declared 1 for multiple Tx streams.
ii) MMSE decoder is used for decoding(calculated in line 115) and (used for combining in line 185 to 191)
iii)ber3 is TXA comparison
iv)ber4 is for TXB comparison