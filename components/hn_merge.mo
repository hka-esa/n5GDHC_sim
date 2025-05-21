within n5GDHC.components;

model hn_merge
  connectors.fluid_p outlet annotation(
    Placement(transformation(extent = {{-50, -10}, {-30, 10}}), iconTransformation(extent = {{-50, -10}, {-30, 10}})));
  connectors.fluid_p inlet_1 annotation(
    Placement(transformation(extent = {{-10, 30}, {10, 50}}), iconTransformation(extent = {{-10, 30}, {10, 50}})));
  connectors.fluid_p inlet_2 annotation(
    Placement(transformation(extent = {{30, -10}, {50, 10}}), iconTransformation(extent = {{30, -10}, {50, 10}})));
equation
  assert(inlet_1.mdot >= (-1e-5), "hn_tube_m I1 negative massflow", level = AssertionLevel.error);
  assert(inlet_2.mdot >= (-1e-5), "hn_tube_m I2 negative massflow", level = AssertionLevel.error);
  assert(outlet.T <= max(inlet_1.T, inlet_2.T) + 1e-10, "hn_tube_m O temperature computation error", level = AssertionLevel.error);
  outlet.mdot = inlet_1.mdot + inlet_2.mdot;
  if inlet_1.mdot < 1e-5 then
    outlet.T = inlet_2.T;
  elseif inlet_2.mdot < 1e-5 then
    outlet.T = inlet_1.T;
  else
    outlet.T = (inlet_1.mdot * inlet_1.T + inlet_2.mdot * inlet_2.T) / (inlet_1.mdot + inlet_2.mdot);
  end if;
  outlet.p = (inlet_1.p * inlet_1.mdot + inlet_2.p * inlet_2.mdot) / (inlet_1.mdot + inlet_2.mdot);
  annotation(
    Icon(coordinateSystem(preserveAspectRatio = false, extent = {{-60, -60}, {60, 60}}), graphics = {Rectangle(lineThickness = 0.5, extent = {{-40, 40}, {40, -40}}), Line( origin = {-21, 0}, rotation = 180,points = {{3, 0}, {-3, 6}, {-3, -6}, {3, 0}}, thickness = 0.5), Line(points = {{-18, 0}, {22, 0}}, thickness = 0.5), Line(origin = {0, 3.04}, rotation = 270, points = {{3, 0}, {-3, 6}, {-3, -6}, {3, 0}}, thickness = 0.5), Line(points = {{0, 24}, {0, 6}}, thickness = 0.5)}),
    Diagram(coordinateSystem(preserveAspectRatio = false, extent = {{-60, -60}, {60, 60}})));
end hn_merge;
