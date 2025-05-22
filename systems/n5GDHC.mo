/*

This file is part of n5GDHC_sim.

Copyright (c) 2025 Manuel Kollmar, Adrian Bürger, Markus Bohlayer, Marco Braun, Moritz Diehl.
Developed at HS Karlsruhe and IMTEK, University of Freiburg.
All rights reserved.

The BSD 3-Clause License

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

within n5GDHC.systems;

model n5GDHC
  // Parameters
  parameter Integer N = 6 "Number of households";
  parameter Integer N_str_1 = 2 "Number of net connections in string 1";
  parameter Integer N_str_2 = 2 "Number of net connections in string 2";
  parameter Integer N_str_3 = 2 "Number of net connections in string 3";
  parameter Integer N_s = 5 "Number of soil layers";
  parameter Integer N_w = 10 "Number of water layers";
  parameter Integer N_hx = 10 "Number of heat exchanger layers";
  parameter Modelica.Units.SI.Length l_n_is = 10.0 "Length of tube part from is";
  parameter Modelica.Units.SI.Radius r_is = 0.08 "Radius network tube";
  parameter Modelica.Units.SI.Radius t_t_is = 0.02 "Thickness tube";
  parameter Modelica.Units.SI.Length l_n = 20 "Length of tube per network connection";
  parameter Modelica.Units.SI.Radius r_n = 0.7 "Radius network tube";
  parameter Modelica.Units.SI.Radius t_t_n = 0.0146 "Thickness tube";
  parameter Modelica.Units.SI.Length l_h = 10 "Length of tube connection -> network";
  parameter Modelica.Units.SI.Radius r_h = 0.02 "Radius connection tube";
  parameter Modelica.Units.SI.Radius t_t_h = 0.004 "Thickness tube";
  parameter Modelica.Units.SI.Volume V_w = 200.0 "Water volume of the storage";
  parameter Modelica.Units.SI.CoefficientOfHeatTransfer alpha_hx_w = 250.0 "Heat transfer coefficient hx - water";
  parameter String path_bypass_inp = Modelica.Utilities.Files.loadResource("modelica://n5GDHC/input/Bypass.txt");
  parameter String path_PI_T_inp = Modelica.Utilities.Files.loadResource("modelica://n5GDHC/input/PI_T.txt"); 
  // Settings
  parameter Modelica.Units.SI.Time t0 = 29376000.0 "Start of simulation in seconds since the beginning of the year";
  parameter Modelica.Units.NonSI.Temperature_degC T_is_b = 2 "Initial bottom temperature of ice storage";
  parameter Modelica.Units.NonSI.Temperature_degC T_is_t = 3 "Initial top temperature of ice storage";
  parameter Modelica.Units.NonSI.Temperature_degC T_n_ret = 2.4 "Initial return temperature of network";
  parameter Modelica.Units.NonSI.Temperature_degC T_n_sup = 2.6 "initial supply temperature of network";
  parameter Modelica.Units.NonSI.Temperature_degC T_n_init = (T_n_ret + T_n_sup) / 2 "Initial fluid temperature";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_0 = T_n_init "Soil layer 1 temperature";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_bd_min = 5.0 "Min soil temperature";
  parameter Modelica.Units.NonSI.Temperature_degC T_s_bd_max = 20.0 "Max soil temperature";
  parameter Modelica.Units.SI.ThermalConductivity lambda_s = 2.0 "Therm. conductivity soil";
  parameter Real ratio_s_w = 0.2 "Volumetric share of water inside soil";
  // PI-controller
  parameter Modelica.Units.NonSI.Temperature_degC T_summer = 19.5 "VL target temperature summer mode";
  //parameter Modelica.Units.NonSI.Temperature_degC T_winter = 1.0 "VL target temperature winter mode";
  parameter Real K_i = 1e-2;
  Real e "Error for PI-controller";
  Real i_pos;
  components.is is(T_w_1_init = T_is_b, T_w_Nw_init = T_is_t, T_hx_init = T_n_init, V_w = V_w, alpha_hx_w = alpha_hx_w, t0 = t0, N_hx = N_hx, N_w = N_w, N_s = N_s, lambda_s = lambda_s, T_s_0 = (T_is_b + T_is_t) / 2, T_s_bd_min = T_s_bd_min, T_s_bd_max = T_s_bd_max, ratio_s_w = ratio_s_w) annotation(
    Placement(transformation(origin = {119, -40.6667}, extent = {{27, -33.3333}, {-27, 26.6667}})));
  components.mvsc mvsc_bypass annotation(
    Placement(transformation(origin = {48.4, -39}, extent = {{8.4, -9}, {-5.6, 9}}, rotation = -180)));
  components.mvm mvm_bypass annotation(
    Placement(transformation(origin = {48.1016, -64}, extent = {{-9.1523, -8}, {6.10156, 8}})));
  Modelica.Blocks.Sources.CombiTimeTable BypassInp(extrapolation = Modelica.Blocks.Types.Extrapolation.HoldLastPoint, fileName = path_bypass_inp, smoothness = Modelica.Blocks.Types.Smoothness.ConstantSegments, tableName = "Bypass", tableOnFile = true) annotation(
    Placement(transformation(origin = {76, 14}, extent = {{10, -10}, {-10, 10}})));
  components.hn_tubes hn_tubes_is(N_s = N_s, l = l_n, r_t = r_is, t_t = t_t_is, t0 = t0, T_f_init = T_n_init, T_s_0 = T_s_0, T_s_bd_min = T_s_bd_min, T_s_bd_max = T_s_bd_max, ratio_s_w = ratio_s_w) annotation(
    Placement(transformation(origin = {-13, -52}, extent = {{32, -24}, {-32, 24}})));
  components.hn_split hn_split_k1 annotation(
    Placement(transformation(origin = {-84, -64}, extent = {{6, -6}, {-6, 6}})));
  components.hn_merge hn_merge_k1 annotation(
    Placement(transformation(origin = {-56, -40}, extent = {{6, -6}, {-6, 6}})));
  components.hn_tubes hn_tubes_k1(N_s = N_s, l = l_n, r_t = r_n, t_t = t_t_n, t0 = t0, T_f_init = T_n_init, T_s_0 = T_s_0, T_s_bd_min = T_s_bd_min, T_s_bd_max = T_s_bd_max, ratio_s_w = ratio_s_w) annotation(
    Placement(transformation(origin = {-70, -2}, extent = {{32, -27}, {-32, 27}}, rotation = -90)));
  components.hn_split hn_split_k3 annotation(
    Placement(transformation(origin = {-84, 40}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  components.hn_merge hn_merge_k3 annotation(
    Placement(transformation(origin = {-56, 40}, extent = {{-6, -6}, {6, 6}}, rotation = 90)));
  constructs.hn_con_h hn_con_h_str_1[N_str_1](each N_s = N_s, each l_n = l_n, each l_h = l_h, each r_n = r_n, each r_h = r_h, each lambda_s = lambda_s, each t0 = t0, each T_init = T_n_init, each t_t_n = t_t_n, each t_t_h = t_t_h, each T_s_0 = T_s_0, each T_s_bd_min = T_s_bd_min, each T_s_bd_max = T_s_bd_max, each ratio_s_w = ratio_s_w);
  constructs.house_param_extr house_param_extr_str_1[N_str_1](idx_r = linspace(1, N_str_1, N_str_1), each T_init = T_n_init);
  constructs.hn_con_h hn_con_h_str_2[N_str_2](each N_s = N_s, each l_n = l_n, each l_h = l_h, each r_n = r_n, each r_h = r_h, each lambda_s = lambda_s, each t0 = t0, each T_init = T_n_init, each t_t_n = t_t_n, each t_t_h = t_t_h, each T_s_0 = T_s_0, each T_s_bd_min = T_s_bd_min, each T_s_bd_max = T_s_bd_max, each ratio_s_w = ratio_s_w);
  constructs.house_param_extr house_param_extr_str_2[N_str_2](idx_r = linspace(N_str_1 + 1, N_str_1 + N_str_2, N_str_2), each T_init = T_n_init);
  constructs.hn_con_h hn_con_h_str_3[N_str_3](each N_s = N_s, each l_n = l_n, each l_h = l_h, each r_n = r_n, each r_h = r_h, each lambda_s = lambda_s, each t0 = t0, each T_init = T_n_init, each t_t_n = t_t_n, each t_t_h = t_t_h, each T_s_0 = T_s_0, each T_s_bd_min = T_s_bd_min, each T_s_bd_max = T_s_bd_max, each ratio_s_w = ratio_s_w);
  constructs.house_param_extr house_param_extr_str_3[N_str_3](idx_r = linspace(N_str_1 + N_str_2 + 1, N_str_1 + N_str_2 + N_str_3, N_str_3), each T_init = T_n_init);
  components.mvsc mvsc_pi annotation(
    Placement(transformation(origin = {76.4, -28}, extent = {{-8.4, -8}, {5.6, 8}})));
  components.mvm mvm_pi annotation(
    Placement(transformation(origin = {77.2, -64}, extent = {{-7.2, -8}, {4.8, 8}})));
  components.pump_c pump_c annotation(
    Placement(transformation(origin = {28, -64}, extent = {{-8, -8}, {8, 8}})));
equation
  connect(hn_split_k1.outlet_2, hn_con_h_str_1[1].inlet_hn_sup);
  connect(hn_con_h_str_1[1].outlet_hn_ret, hn_merge_k1.inlet_2);
  for k in 1:N_str_1 - 1 loop
    connect(hn_con_h_str_1[k].outlet_h_sup, house_param_extr_str_1[k].inlet_hn);
    connect(house_param_extr_str_1[k].outlet_hn, hn_con_h_str_1[k].inlet_h_ret);
    connect(hn_con_h_str_1[k].outlet_hn_sup, hn_con_h_str_1[k + 1].inlet_hn_sup);
    connect(hn_con_h_str_1[k + 1].outlet_hn_ret, hn_con_h_str_1[k].inlet_hn_ret);
  end for;
  connect(hn_con_h_str_1[N_str_1].outlet_h_sup, house_param_extr_str_1[N_str_1].inlet_hn);
  connect(house_param_extr_str_1[N_str_1].outlet_hn, hn_con_h_str_1[N_str_1].inlet_h_ret);
  hn_con_h_str_1[N_str_1].outlet_hn_sup.mdot = 0.0;
  hn_con_h_str_1[N_str_1].inlet_hn_ret.mdot = 0.0;
  hn_con_h_str_1[N_str_1].inlet_hn_ret.T = 0.0;
  hn_con_h_str_1[N_str_1].inlet_hn_ret.p = 0.0;
  connect(hn_split_k3.outlet_1, hn_con_h_str_2[1].inlet_hn_sup);
  connect(hn_con_h_str_2[1].outlet_hn_ret, hn_merge_k3.inlet_1);
  for k in 1:N_str_2 - 1 loop
    connect(hn_con_h_str_2[k].outlet_h_sup, house_param_extr_str_2[k].inlet_hn);
    connect(house_param_extr_str_2[k].outlet_hn, hn_con_h_str_2[k].inlet_h_ret);
    connect(hn_con_h_str_2[k].outlet_hn_sup, hn_con_h_str_2[k + 1].inlet_hn_sup);
    connect(hn_con_h_str_2[k + 1].outlet_hn_ret, hn_con_h_str_2[k].inlet_hn_ret);
  end for;
  connect(hn_con_h_str_2[N_str_2].outlet_h_sup, house_param_extr_str_2[N_str_2].inlet_hn);
  connect(house_param_extr_str_2[N_str_2].outlet_hn, hn_con_h_str_2[N_str_2].inlet_h_ret);
  hn_con_h_str_2[N_str_2].outlet_hn_sup.mdot = 0.0;
  hn_con_h_str_2[N_str_2].inlet_hn_ret.mdot = 0.0;
  hn_con_h_str_2[N_str_2].inlet_hn_ret.T = 0.0;
  hn_con_h_str_2[N_str_2].inlet_hn_ret.p = 0.0;
  connect(hn_split_k3.outlet_2, hn_con_h_str_3[1].inlet_hn_sup);
  connect(hn_con_h_str_3[1].outlet_hn_ret, hn_merge_k3.inlet_2);
  for k in 1:N_str_3 - 1 loop
    connect(hn_con_h_str_3[k].outlet_h_sup, house_param_extr_str_3[k].inlet_hn);
    connect(house_param_extr_str_3[k].outlet_hn, hn_con_h_str_3[k].inlet_h_ret);
    connect(hn_con_h_str_3[k].outlet_hn_sup, hn_con_h_str_3[k + 1].inlet_hn_sup);
    connect(hn_con_h_str_3[k + 1].outlet_hn_ret, hn_con_h_str_3[k].inlet_hn_ret);
  end for;
  connect(hn_con_h_str_3[N_str_3].outlet_h_sup, house_param_extr_str_3[N_str_3].inlet_hn);
  connect(house_param_extr_str_3[N_str_3].outlet_hn, hn_con_h_str_3[N_str_3].inlet_h_ret);
  hn_con_h_str_3[N_str_3].outlet_hn_sup.mdot = 0.0;
  hn_con_h_str_3[N_str_3].inlet_hn_ret.mdot = 0.0;
  hn_con_h_str_3[N_str_3].inlet_hn_ret.T = 0.0;
  hn_con_h_str_3[N_str_3].inlet_hn_ret.p = 0.0;
// PI
  e = if time < 6249600 then T_summer - mvsc_bypass.inlet.T else mvsc_bypass.inlet.T - T_summer;
  der(i_pos) = if (i_pos >= 1 and e > 0) or (i_pos <= 0 and e < 0) then 0 else K_i*e;
  mvsc_pi.pos = i_pos;
//
  connect(hn_tubes_k1.outlet_hn_return, hn_merge_k1.inlet_1) annotation(
    Line(points = {{-56.5, -23.3333}, {-56.5, -35.3333}}, color = {28, 108, 200}));
  connect(hn_split_k1.outlet_1, hn_tubes_k1.inlet_hn_supply) annotation(
    Line(points = {{-84, -60}, {-84, -24}}, color = {28, 108, 200}));
  connect(hn_tubes_k1.outlet_hn_supply, hn_split_k3.inlet) annotation(
    Line(points = {{-84, 20}, {-84, 36}}, color = {28, 108, 200}));
  connect(hn_merge_k3.outlet, hn_tubes_k1.inlet_hn_return) annotation(
    Line(points = {{-56, 36}, {-56, 20}}, color = {28, 108, 200}));
  connect(hn_merge_k1.outlet, hn_tubes_is.inlet_hn_return) annotation(
    Line(points = {{-52, -40}, {-34, -40}}, color = {28, 108, 200}));
  connect(hn_tubes_is.outlet_hn_supply, hn_split_k1.inlet) annotation(
    Line(points = {{-34, -64}, {-80, -64}}, color = {28, 108, 200}));
  connect(hn_tubes_is.outlet_hn_return, mvsc_bypass.inlet) annotation(
    Line(points = {{8, -40}, {43, -40}, {43, -39}}, color = {28, 108, 200}));
  connect(mvsc_bypass.outlet_2, mvsc_pi.inlet) annotation(
    Line(points = {{48.4, -33}, {48.4, -29}, {70.4, -29}}, color = {28, 108, 200}));
  connect(mvsc_bypass.outlet_1, mvm_bypass.inlet_1) annotation(
    Line(points = {{48.4, -45}, {48.4, -60}}, color = {28, 108, 200}));
  connect(mvsc_pi.outlet_2, mvm_pi.inlet_1) annotation(
    Line(points = {{76.4, -33.3333}, {78.4, -33.3333}, {78.4, -57.3333}}, color = {28, 108, 200}));
  connect(mvm_pi.outlet, mvm_bypass.inlet_2) annotation(
    Line(points = {{72.4, -64}, {60.4, -64}, {60.4, -78}, {48.4, -78}, {48.4, -70}}, color = {28, 108, 200}));
  connect(mvsc_pi.outlet_1, is.inlet_hx) annotation(
    Line(points = {{76.4, -22.6667}, {76.4, -14.6667}, {88.4, -14.6667}, {88.4, -34.6667}, {98.4, -34.6667}}, color = {28, 108, 200}));
  connect(is.outlet_hx, mvm_pi.inlet_2) annotation(
    Line(points = {{98.75, -61.3333}, {88.75, -61.3333}, {88.75, -75.3333}, {78.75, -75.3333}, {78.75, -69.3333}}, color = {28, 108, 200}));
  connect(BypassInp.y[1], mvsc_bypass.pos) annotation(
    Line(points = {{66, 14}, {56, 14}, {56, -38}, {50, -38}}, color = {0, 0, 127}));
  connect(mvm_bypass.outlet, pump_c.inlet) annotation(
    Line(points = {{42, -64}, {34, -64}}, color = {28, 108, 200}));
  connect(pump_c.outlet, hn_tubes_is.inlet_hn_supply) annotation(
    Line(points = {{22, -64}, {8, -64}}, color = {28, 108, 200}));
  annotation(
    Diagram);
end n5GDHC;
