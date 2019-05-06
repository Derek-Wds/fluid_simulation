/// <reference path="./p5.global-mode.d.ts" />

/*
    Function of diffuse
    - b : int
    - x : float[]
    - x0 : float[]
    - diff : float
    - dt : flaot
*/
function diffuse(b, x, x0, diff, dt) {
    let a = dt * diff * (N - 2) * (N - 2);
    lin_solve(b, x, x0, a, 1 + 6 * a);
}


