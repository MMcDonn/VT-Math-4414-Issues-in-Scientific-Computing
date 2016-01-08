function w = gaussian_quadr_weights(a, b, n)
bw = gaussian_quadr_weights_ref(n);
w = ((b-a)/2)*bw;
return;