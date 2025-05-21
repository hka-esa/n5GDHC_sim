within n5GDHC.components;

model pump_c
  parameter Modelica.Units.SI.Efficiency eta = 0.912 "Motor efficiency at full power";
  parameter Modelica.Units.SI.Pressure p_ref = 3.4/1e-5 "Reference pressure for supply line";
  parameter Modelica.Units.SI.Density rho_f = 1.0e3 "Density of network medium";
  Modelica.Units.SI.Pressure p_diff "Pressure increase";
  Modelica.Units.SI.Power P_el "Electric consumption of pump";
  Modelica.Units.SI.Energy E_el(start=0, fixed=true) "Energy consumed";
  connectors.fluid_p outlet annotation(
    Placement(transformation(extent = {{-70, -10}, {-50, 10}})));
  connectors.fluid_p inlet annotation(
    Placement(transformation(extent = {{50, -10}, {70, 10}})));
equation
  assert(abs(outlet.mdot - inlet.mdot) <= 1e-5, "Pump mass balance violated.", level = AssertionLevel.error);
  outlet.mdot = inlet.mdot;
  outlet.T = inlet.T;
  outlet.p = p_ref;
  p_diff = p_ref - inlet.p;  
  P_el = p_diff * inlet.mdot / (rho_f * eta);
  der(E_el) = P_el;  
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -80}, {80, 80}}), graphics = {Ellipse(lineThickness = 0.5, extent = {{-60, 60}, {60, -60}}), Line(points = {{-56, 20}, {20, 56}}, thickness = 0.5), Line(points = {{-56, -20}, {20, -56}}, thickness = 0.5)}),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-80, -80}, {80, 80}})));
end pump_c;
