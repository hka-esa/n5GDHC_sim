within n5GDHC.components;

model mvm
  connectors.fluid_p outlet annotation(
    Placement(transformation(extent = {{-50, -10}, {-30, 10}}), iconTransformation(extent = {{-50, -10}, {-30, 10}})));
  connectors.fluid_p inlet_1 annotation(
    Placement(transformation(extent = {{-10, 30}, {10, 50}}), iconTransformation(extent = {{-10, 30}, {10, 50}})));
  connectors.fluid_p inlet_2 annotation(
    Placement(transformation(extent = {{-10, -50}, {10, -30}}), iconTransformation(extent = {{-10, -50}, {10, -30}})));
equation
  assert(inlet_1.mdot >= (-1e-5), "MVC I1 negative massflow", level = AssertionLevel.error);
  assert(inlet_2.mdot >= (-1e-5), "MVC I2 negative massflow", level = AssertionLevel.error);
  assert(outlet.T - 1e-5 <= max(inlet_1.T, inlet_2.T), "MV O temperature computation error", level = AssertionLevel.error);
  outlet.mdot = inlet_1.mdot + inlet_2.mdot;
  if inlet_1.mdot < 1e-5 then
    outlet.T = inlet_2.T;
    outlet.p = inlet_2.p;
  elseif inlet_2.mdot < 1e-5 then
    outlet.T = inlet_1.T;
    outlet.p = inlet_1.p;
  else
    outlet.T = (inlet_1.mdot * inlet_1.T + inlet_2.mdot * inlet_2.T) / (inlet_1.mdot + inlet_2.mdot);
    outlet.p = (inlet_1.p * inlet_1.mdot + inlet_2.p * inlet_2.mdot) / (inlet_1.mdot + inlet_2.mdot);
  end if;
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-60, -60}, {40, 60}}), graphics = {Line(points = {{-20, 40}, {20, -40}, {-20, -40}, {20, 40}, {-20, 40}}, thickness = 0.5), Line(points = {{0, 0}, {-40, 20}, {-40, -20}, {0, 0}}, thickness = 0.5), Line(points = {{-18, 6}, {-18, -6}, {-26, 0}, {-18, 6}}, thickness = 0.5)}),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-60, -60}, {40, 60}})));
end mvm;
