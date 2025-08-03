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

within n5GDHC.components;

/* 
The parameter data is taken from literature, public sources or data sheets in order to reflect the specific conditions at the site.
Sources:
[1] Kuchling, H.: Taschenbuch der Physik, Hanser München, 22, 2022.
*/

model hx_param_qdot
  // States and initials
  parameter Modelica.Units.NonSI.Temperature_degC T_init = 10.0 "Initial temperature";
  Modelica.Units.NonSI.Temperature_degC T_B(start = T_init, fixed = true) "Temperatures side B";
  // Material properties
  parameter Modelica.Units.SI.SpecificHeatCapacity c_B = 3.8e3 "Specific heat capacity medium side B";
  parameter Modelica.Units.SI.Density rho_B = 1.0e3 "Density medium side B [1]";
  // Geometrical propteries
  parameter Modelica.Units.SI.Volume V = 0.005 "Volume per side";
  // Calculated variables
  Modelica.Units.SI.Mass m_B "Mass medium side B";
  // Monitoring
  Modelica.Units.SI.HeatFlowRate Qdot_hx;
  Modelica.Units.SI.Heat Q_hx_inp(start = 0, fixed = true);
  Modelica.Units.SI.Heat Q_hx(start = 0, fixed = true);
  connectors.fluid_p inlet_B annotation(
        Placement(transformation(extent = {{20, 30}, {40, 50}}), iconTransformation(extent = {{-40, -50}, {-20, -30}})));
  connectors.fluid_p outlet_B annotation(
        Placement(transformation(extent = {{30, 30}, {50, 50}}), iconTransformation(extent = {{20, -50}, {40, -30}})));
  Modelica.Blocks.Interfaces.RealInput Qdot_ext annotation(
    Placement(visible = true, transformation(origin = {-20, 24}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-3, 37}, extent = {{-13, -13}, {13, 13}}, rotation = -90)));
equation
  assert(inlet_B.mdot >= (-1e-5), "HX IB negative massflow", level = AssertionLevel.error);
  m_B = V * rho_B;
  m_B * c_B * der(T_B) = inlet_B.mdot * c_B * (inlet_B.T - T_B) + Qdot_ext;
  outlet_B.mdot = inlet_B.mdot;
  outlet_B.T = T_B;
// Monitoring
  Qdot_hx = inlet_B.mdot * c_B * (outlet_B.T - inlet_B.T);
  der(Q_hx_inp) = Qdot_ext;
  der(Q_hx) = Qdot_hx;
// pressure
  inlet_B.p = outlet_B.p;
  annotation(
        Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -60}, {80, 60}})),
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -60}, {80, 60}}), graphics = {Rectangle(fillColor = {255, 255, 255}, fillPattern = FillPattern.Solid, lineThickness = 0.5, extent = {{-60, 40}, {60, -40}}), Line(points = {{-60, 40}, {60, -40}}, thickness = 0.5), Text( extent = {{26, 12}, {58, -10}}, textString = "A"), Text( extent = {{-60, 12}, {-28, -10}}, textString = "B")}),
    experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002));
end hx_param_qdot;
