# nonAdiabatic
Check nonadiabatic spin flip loss in a pulsed decelerator

I want to study non-adiabatic transitions in a decelerator. I could do
this by getting the E-field magnitude and direction along the trajectory
of various molecules of interest. I think it will be reasonable to just
assume the trajectories are parallel to the central axis, since their
velocity is so dominant in this direction initially.

I cannot simply use the data stored in the slowANDtrap simulation, which
only has the internal energy of the molecule. I could invert this to get
E-field, but I wouldn't know the field direction, which is definitely
relevant for f-3/2 to f-1/2 transitions, which are only coupled by a
rotating E-field.

In fact, the only COMSOL generated data that does have true field
direction and not only the E-B angle is a part of the PinTrapLZ
simulations, which are not close enough to the actual decelerator
geometry.

So it's back to COMSOL. And what will I export? I'd prefer to have a grid
of lines parallel to the axis, which entails a 3D array. It's not enough
to just start along the axis, because this might have an unusual
magnitude of rotation due to the symmetry.
