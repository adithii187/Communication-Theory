QPSK Type 2 Communication Model Simulation in MATLAB
Introduction

This project involves simulating a complete communication model using MATLAB. The primary goal is to recreate a given audio file after processing it through various stages of the communication model. The project covers the following key blocks: A/D converter, Encoder, Line Coding, Modulation, Channel, and Demodulation. The specific encoding scheme implemented is Quadrature Phase-Shift Keying (QPSK) Type 2.
Procedure
1. A/D Converter

    Objective: Convert the provided .wav audio file into a bitstream.
    Method: Utilize MATLAB functions to read the audio file and convert it to a binary bitstream.
    Result: Bitstream representing the audio signal.

2. Encoder

    Objective: Convert the bitstream into a sequence of symbols.
    Encoding Scheme: QPSK Type 2.
        Constellation Mapping: Based on the constellation diagram for QPSK Type 2 (as shown in Figure 2(b) of the provided documentation).
    Method: Map pairs of bits to corresponding QPSK symbols.

3. Line Coding

    Objective: Generate the line-coded signal.
    Equation: x3(t)=∑kakp(t−kTb)x3​(t)=∑k​ak​p(t−kTb​), where p(t)p(t) is a pulse of duration TbTb​.
    Pulse Shapes:
        Rectangular pulse.
        Raised cosine pulse.
    Result: Line-coded signal for both pulse shapes.

4. Modulation

    Objective: Modulate the line-coded signal.
    Modulation Scheme: Phase modulation for QPSK.
    Implementation: Use MATLAB functions to modulate the signal using QPSK.
    Carrier Frequency: 100 Hz for visualization; 1 MHz for final communication model.
    Result: Modulated signal ready for transmission through the channel.

5. Channel and Demodulation

    Objective: Transmit the signal through a channel and demodulate it.
    Channels:
        Memoryless AWGN channel: r(t)=s(t)+n(t)r(t)=s(t)+n(t).
        AWGN channel with memory: r(t)=h(t)∗s(t)+n(t)r(t)=h(t)∗s(t)+n(t), where h(t)=aδ(t)+(1−a)δ(t−bTb)h(t)=aδ(t)+(1−a)δ(t−bTb​).
    Demodulation: Process the received signal to recover the original bitstream.
    Result: Demodulated signal, which is then decoded back to the original audio.

Outputs

    Plot the output of every block:
        x1(t)x1​(t): Output of the A/D converter.
        x2(t)x2​(t): Output of the encoder.
        x3(t)x3​(t): Output of the line coding block.
        x4(t)x4​(t): Output of the modulation block.
        y1(t)y1​(t): Output after passing through the memoryless AWGN channel.
        y2(t)y2​(t): Output after demodulation from the memoryless AWGN channel.
        y3(t)y3​(t): Output after passing through the AWGN channel with memory.
        y4(t)y4​(t): Output after demodulation from the AWGN channel with memory.

    Plot the Power Spectral Density (PSD) of:
        x3(t)x3​(t)
        x4(t)x4​(t)
        y3(t)y3​(t)
        y4(t)y4​(t)

    Probability of Error PePe​ vs SNR:
        For both channel realizations (memoryless AWGN and AWGN with memory).
        Provide comparative analysis.

    Constellation Diagrams:
        Input and output constellations for both channel realizations.
        Observe and analyze the differences.

Conclusion

This project demonstrates the end-to-end simulation of a communication model using QPSK Type 2 encoding. The various stages from A/D conversion to channel transmission and demodulation have been successfully implemented and analyzed. The results include plots and analysis of each block, error probabilities, and constellation diagrams, providing a comprehensive understanding of the communication process.
Files Included

    main.m: MATLAB script containing the main simulation code.
    functions.m: MATLAB functions used for various blocks.
    audio.wav: Original audio file to be processed.
    report.pdf: Detailed report including plots, calculations, and analysis.

Usage

    Run main.m to execute the entire communication model simulation.
    Refer to report.pdf for detailed explanations and analysis of the results.

References

    MATLAB Documentation for signal processing functions.
    Course materials on digital communication systems and modulation techniques.
