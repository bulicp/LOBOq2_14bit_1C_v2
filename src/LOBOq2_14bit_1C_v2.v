// Qx = 8, Qy = 12

module LOBOq2_14bit_1C_v2(
    input [15:0] x,
    input [15:0] y,
    output [31:0] p
    );

    // Generate pp
    
    // Radix 4 pp
    
    wire [16:0] pp_rad4_0;
    wire sign_factor;
    
    rad4_gen12 gen_rad4_pp( 
        .x1(x[15:13]), .y(y),
        .pp_rad4_0(pp_rad4_0),
        .sign_factor(sign_factor));
        
    // Radix 1024 pp 
   
    wire [16:0] PP_oa;
    wire [28:0] PP_ob;
    
    rad4096_gen_v1 gen_rad4096_pp( 
      .x2(x[13:0]),.y(y),
      .PP_oa(PP_oa),
      .PP_ob(PP_ob));
   
    // Tree
    pp_tree_red_v2_12bit wallace_tree(
        .pp_rad4_0(pp_rad4_0),
        .PP_oa(PP_oa),
        .PP_ob(PP_ob),
        .sign_factor(sign_factor),
        .p(p));

endmodule

module rad4_gen12( x1,y,pp_rad4_0,sign_factor);
// inputs
// y multiplicand
// x multipland 
// P1,P2,P3 partial products
input [15:0] y;
input [2:0] x1;
// output
output [16:0] pp_rad4_0;
output  sign_factor;

wire one,two,sign;


code code0(one,two,sign,x1[2],x1[1], x1[0]);



wire [16:0] tmp1_pp; 
assign tmp1_pp = {y[15],y}; // This variable is introduced because pp has 17 bits

// first pp generation 
wire [17:0] out1;
assign out1[0] = sign;

genvar i;
generate
	for ( i = 0; i < 17; i = i+1 )
		begin : pp_rad4_first 
		product pp_pr(tmp1_pp[i],out1[i],one,two,sign,pp_rad4_0[i],out1[i+1]);
		end
endgenerate


sgn_gen sgn_x(one,two,sign,sign_factor);
endmodule

module code(one,two,sign,y2,y1,y0);
	input y2,y1,y0;
	output one,two,sign;
	wire [1:0]k;
	xor x1(one,y0,y1);
	xor x2(k[1],y2,y1);
	not n1(k[0],one);
	and a1(two,k[0],k[1]);
	assign sign=y2;
endmodule

module product(x1,x0,one,two,sign,p,i);
	input x1,x0,sign,one,two;
	output p,i;
	wire [1:0] k;
	xor xo1(i,x1,sign);
	and a1(k[1],i,one);
	and a0(k[0],x0,two);
	or o1(p,k[1],k[0]);
endmodule

module sgn_gen(one,two,sign,sign_factor);
	input sign,one,two;
	output sign_factor;
	wire k;
	or o1(k,one,two);
	and a1(sign_factor,sign,k);
endmodule

module rad4096_gen_v1(x2,y,PP_oa,PP_ob);

// inputs
// y multiplicand
// x multipland 
// P0a -> x,P0b -> partial products
input [15:0] y;
input [13:0] x2;
// output
output [16:0] PP_oa;
output [28:0] PP_ob;


// LOD detection X
// First calculate abs(x2)
// then apply LOD

wire [13:0] abs_x2;
wire [6:0] k_x2;
wire [2:0] k_x2_enc;
wire zero_x2;


sec_complement_w14 x0_abs_gen
		(
            .data_in(x2), 
            .sign(x2[13]),
            .data_out(abs_x2) 
        );

LOD14_quant lod_x2(abs_x2,zero_x2,k_x2);

PriorityEncoder_7 PEx(k_x2,k_x2_enc);

// LOD detection Y
wire [15:0] y_abs;
wire [3:0] k_y;
wire [1:0] k_y_enc;
wire zero_y;

sec_complement_w16 y_abs_gen 
		(
            .data_in(y), 
            .sign(y[15]),
            .data_out(y_abs) 
        );	


LOD16_quant lod16_y(y_abs[15:12],zero_y,k_y);
PriorityEncoder_4 PEy(k_y,k_y_enc);


// Reset bit on k_x position of X and invert it i
wire [13:0] x20; // value after bit reset


assign x20[13:8] =  abs_x2[13:8] & (~k_x2[6:1]);
assign x20[7:1] = abs_x2[7:1];
assign x20[0] = abs_x2[0] & (~k_x2[0]);

// generation of P0a and P0b

// generation of P0a and P0b
wire prod_sign;
wire [13:0] tmp_out0;
wire [15:0] tmp_out1;

assign prod_sign = x2[13] ^ y[15];


wire [13:0] x20_signed;
wire [15:0] y_signed;

sec_complement_w14 x20_sign_gen 
		(
            .data_in(x20), 
            .sign(prod_sign),
            .data_out(tmp_out0) 
        );

		  
sec_complement_w16 y_sign_gen 
		(
            .data_in(y_abs), 
            .sign(prod_sign),
            .data_out(tmp_out1) 
        );

assign y_signed = (!zero_x2) ? tmp_out1 : 16'b0;			
assign x20_signed = (!zero_y) ? tmp_out0 : 14'b0;


Barrel27L_14_quant gen_P0b(
			  .shift_i(k_x2_enc),
             .data_i(y_signed),
             .data_o(PP_ob));
						  
Barrel15L_16 gen_P0a(
            .data_i(x20_signed),
            .shift_i(k_y_enc),
            .data_o(PP_oa));
			
endmodule 

module LOD16_quant(
    input [3:0] data_i,
    output zero_o,
    output [3:0] data_o
    );
	 
	 wire [3:0] z;
	 
	 //*****************************************
	 // Zero detection logic:
	 //*****************************************
	 assign zdet = data_i[3] | data_i[2] | data_i[1] | data_i[0];
	 assign zero_o = ~ zdet;
		 
		 
	 //*****************************************
	 // LODs:
	 //*****************************************
	 LOD4 lod4_1 (
		.data_i(data_i), 
		.data_o(data_o)
	 );
	 
	 
	 //*****************************************
	 // Multiplexers :
	 //*****************************************

endmodule


module LOD14_quant(
    input [13:0] data_i,
    output zero_o,
    output [6:0] data_o
    );
	
	 wire [5:0] tmp_in;
	 wire [5:0] z;
	 wire [1:0] zdet;
	 wire [1:0] select;
	 //*****************************************
     // Quantization 
     //*****************************************
	 assign tmp_in[5:4] = data_i[12:11];
	 assign tmp_in[3:1] = data_i[10:8];
	 assign tmp_in[0] =  |data_i[7:0];
	 //*****************************************
	 // Zero detection logic:
	 //*****************************************
	 assign zdet[0] = tmp_in[3] | tmp_in[2] | tmp_in[1] | tmp_in[0];
	 assign zdet[1] = tmp_in[5] | tmp_in[4];
	 assign zero_o = ~(data_i[13] |  zdet[0]  | zdet[1] );
		 
		 
	 //*****************************************
	 // LODs:
	 //*****************************************
	 LOD4 lod4_1 (
		.data_i(tmp_in[3:0]), 
		.data_o(z[3:0])
	 );
	 
	  LOD2 lod2_2 (
            .data_i(tmp_in[5:4]), 
            .data_o(z[5:4])
         );
	 LOD2 lod2_3 (
                .data_i(zdet), 
                .data_o(select)
             );
	 
	 //*****************************************
	 // Multiplexers :
	 //*****************************************	
	 assign data_o[6] = data_i[13];
	 
	 Muxes2in1Array2 Inst_MUX214_1 (
        .data_i(z[5:4]), 
        .select_i(select[1]), 
        .data_o(data_o[5:4])
    );

	 Muxes2in1Array4 Inst_MUX214_0 (
		.data_i(z[3:0]), 
		.select_i(select[0]), 
		.data_o(data_o[3:0])
    );

endmodule


module LOD3(
    input [2:0] data_i,
    output [2:0] data_o
    );
	 
	 
	 wire mux0;
	 wire mux1;
	
	 
	 // multiplexers:
	 assign mux1 = (data_i[2]==1) ? 1'b0 : 1'b1;
	 assign mux0 = (data_i[1]==1) ? 1'b0 : mux1;
	 
	 //gates and IO assignments:
	 assign data_o[2] = data_i[2];
	 assign data_o[1] =(mux1 & data_i[1]);
	 assign data_o[0] =(mux0 & data_i[0]);
	 
endmodule

module LOD4(
    input [3:0] data_i,
    output [3:0] data_o
    );
	 
	 
	 wire mux0;
	 wire mux1;
	 wire mux2;
	 
	 // multiplexers:
	 assign mux2 = (data_i[3]==1) ? 1'b0 : 1'b1;
	 assign mux1 = (data_i[2]==1) ? 1'b0 : mux2;
	 assign mux0 = (data_i[1]==1) ? 1'b0 : mux1;
	 
	 //gates and IO assignments:
	 assign data_o[3] = data_i[3];
	 assign data_o[2] =(mux2 & data_i[2]);
	 assign data_o[1] =(mux1 & data_i[1]);
	 assign data_o[0] =(mux0 & data_i[0]);
	 

endmodule

module LOD2(
    input [1:0] data_i,
    output [1:0] data_o
    );
	 
	 
	 //gates and IO assignments:
	 assign data_o[1] = data_i[1];
	 assign data_o[0] =(~data_i[1] & data_i[0]);
	 

endmodule

module Muxes2in1Array4(
    input [3:0] data_i,
    input select_i,
    output [3:0] data_o
    );

	assign data_o[3] = select_i ? data_i[3] : 1'b0;
	assign data_o[2] = select_i ? data_i[2] : 1'b0;
	assign data_o[1] = select_i ? data_i[1] : 1'b0;
	assign data_o[0] = select_i ? data_i[0] : 1'b0;
	
endmodule

module Muxes2in1Array2(
    input [1:0] data_i,
    input select_i,
    output [1:0] data_o
    );
    
	assign data_o[1] = select_i ? data_i[1] : 1'b0;
	assign data_o[0] = select_i ? data_i[0] : 1'b0;
	
endmodule

module PriorityEncoder_4(
    input [3:0] data_i,
    output reg [1:0] code_o
    );

	always @*
		case (data_i)
	     4'b0001 : code_o = 3'b00;
         4'b0010 : code_o = 3'b01;
         4'b0100 : code_o = 3'b10;
		 default  : code_o = 3'b11;
		endcase		
endmodule

module PriorityEncoder_7(
    input [6:0] data_i,
    output reg [2:0] code_o
    );

	  always @*
		case (data_i)
	     7'b0000001 : code_o = 3'b000;
         7'b0000010 : code_o = 3'b001;
         7'b0000100 : code_o = 3'b010;
         7'b0001000 : code_o = 3'b011;
         7'b0010000 : code_o = 3'b100;
         7'b0100000 : code_o = 3'b101;
		 default  : code_o = 3'b110;
		endcase		
	endmodule



module Barrel27L_14_quant(
    input [15:0] data_i,
    input [2:0] shift_i,
    output reg [28:0] data_o
    );
	
	wire [28:0] tmp;
   assign tmp = {{13{data_i[15]}},data_i};
  // assign tmp = {data_i};
  always @*
     case (shift_i)
        3'b000: data_o = tmp;
        3'b001: data_o = tmp << 8;
        3'b010: data_o = tmp << 9;
	    3'b011: data_o = tmp << 10;
	    3'b100: data_o = tmp << 11;		
        default: data_o = tmp << 12;
     endcase
endmodule

module Barrel15L_16(
    input [13:0] data_i,
    input [1:0] shift_i,
    output reg [16:0] data_o
    );
   wire [16:0] tmp;
   assign tmp = {{3{data_i[13]}},data_i};
  //assign tmp = {data_i};
   always @*
      case (shift_i)
         2'b00: data_o = tmp;
         2'b01: data_o = tmp << 1;
         2'b10: data_o = tmp << 2;
         default: data_o = tmp << 3;
      endcase


endmodule


module pp_tree_red_v2_12bit(pp_rad4_0,PP_oa,PP_ob,sign_factor,p);
// inputs
// pp_rad4_x - Rad4
// PP_oa, PP_ob - Rad1024
input [16:0] pp_rad4_0;
input [16:0] PP_oa;
input [28:0] PP_ob;
input  sign_factor;
// output
// product p
output [31:0] p;

// generate negative MSBs
wire [2:0] E_MSB;
not n1(E_MSB[0],PP_ob[28]);
not n2(E_MSB[1],PP_oa[16]);
not n3(E_MSB[2],pp_rad4_0[16]);


// Reduction 

// First reduction

// first group
wire [15:0] sum00_FA;
wire [15:0] carry00_FA;


wire [15:0] tmp001_FA;
wire [15:0] tmp002_FA;
wire [15:0] tmp003_FA;

assign tmp001_FA = {E_MSB[0],PP_ob[28:14]};
assign tmp002_FA = {E_MSB[1],PP_oa[16:2]};
assign tmp003_FA = {pp_rad4_0[15:0]};


genvar i001;
generate
	for (i001 = 0; i001 < 16; i001 = i001 + 1)
		begin : pp_fad00
		FAd pp_fad(tmp001_FA[i001],tmp002_FA[i001], tmp003_FA[i001], carry00_FA[i001],sum00_FA[i001]);
		end
endgenerate

wire sum00_HA;
wire carry00_HA;

assign sum00_HA = E_MSB[2];
assign carry00_HA = pp_rad4_0[16];

// CLA adder 

wire [31:0] tmp_ADD1;
wire [31:0] tmp_ADD2; 

assign tmp_ADD1 = {E_MSB[2],sum00_HA, sum00_FA,PP_ob[13:0]};
assign tmp_ADD2 = {carry00_HA, carry00_FA,sign_factor,PP_oa[1:0],12'b0};

// Final product

assign p = tmp_ADD1 + tmp_ADD2;

endmodule 

//adders design
module HAd(a,b,c,s);
	input a,b;
	output c,s;
	xor x1(s,a,b);
	and a1(c,a,b);
endmodule

module FAd(a,b,c,cy,sm);
	input a,b,c;
	output cy,sm;
	wire x,y,z;
	xor x1(x,a,b);
	xor x2(sm,x,c);
	and a1(y,a,b);
	and a2(z,x,c);
	or o1(cy,y,z);
endmodule

module FA(a,b,c,cy,sm);
	input a,b,c;
	output cy,sm;
	wire x,y,z;
	xor x1(x,a,b);
	xor x2(sm,x,c);
	and a1(y,a,b);
	and a2(z,x,c);
	or o1(cy,y,z);
endmodule


module sec_complement_w14
  (
   input [13:0] data_in,
   input sign,
   output [13:0] data_out
   );
     

 
  // Create the HA Adders
  genvar  ii;
  generate
    for (ii=0; ii<14; ii=ii+1) 
      begin: pc
         //assign data_out[ii] = data_in[ii] ^ (sign & w_C[ii-1]);
           assign data_out[ii] = data_in[ii] ^ (sign);

      end
  endgenerate
 
  // Create the Generate (G) Terms:  Gi=Ai*Bi
  // Create the Propagate Terms: Pi=Ai+Bi
  // Create the Carry Terms:

 
endmodule // carry_lookahead_adder

module sec_complement_w16
  (
   input [15:0] data_in,
   input sign,
   output [15:0] data_out
   );
     

 
  // Create the HA Adders
  genvar  ii;
  generate
    for (ii=0; ii<16; ii=ii+1) 
      begin: pc
        // assign data_out[ii] = data_in[ii] ^ (sign & w_C[ii-1]);
	   assign data_out[ii] = data_in[ii] ^ (sign);

      end
  endgenerate
 
  // Create the Generate (G) Terms:  Gi=Ai*Bi
  // Create the Propagate Terms: Pi=Ai+Bi
  // Create the Carry Terms:

endmodule // carry_lookahead_adder

