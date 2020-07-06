
`include "timescale.v"

module mult_top(
    // Clock and reset
    input clk,
    // Input X has a form x_i%d where %d denotes the bit number
    input  [15:0] x,
    // Input Y has a form y_i%d where %d denotes the bit number
    input  [15:0] y,
    // Output P has a form p_out%d where %d denotes the bit number
    output reg [31:0] p_out
    );

    
    // Now we have X_vec and Y_vec signal 
    // Then we do processing with these signals and store the 
    // intermidiate result in P_vec
    // For example purposes X_vec and Y_vec are concanated and stored in P_vec
    wire [31:0] P_vec;
    reg [15:0] X_vec;
    reg [15:0] Y_vec;
    
    LOBOq2_14bit_1C_v2 mult(X_vec,Y_vec,P_vec); 

 
    always @(posedge clk) 
    begin
        p_out = P_vec;
        X_vec = x;
        Y_vec = y;
    end

endmodule 

