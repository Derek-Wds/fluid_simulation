/// <reference path="./p5.global-mode.d.ts" />

let N = 256;
let iter = 16;
let SCALE = 4;
let t = 0;

// function to use 1D array and fake the extra two dimensions --> 3D
function IX(x, y, z) {
    return x + y*N + z*N*N;
}

// Fluid cube class
class Fluid {
    constructor(dt, diffusion, viscosity) {
        this.size = N;
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;

        this.s = new Array(N*N*N);
        this.density = new Array(N*N*N);

        this.Vx = new Array(N*N*N);
        this.Vy = new Array(N*N*N);
        this.Vz = new Array(N*N*N);

        this.Vx0 = new Array(N*N*N);
        this.Vy0 = new Array(N*N*N);
        this.Vz0 = new Array(N*N*N);
    }

    // step method
    step() {
        let N = this.size;
        let visc = this.visc;
        let diff = this.diff;
        let dt = this.dt;
        let Vx = this.Vx;
        let Vy = this.Vy;
        let Vz = this.Vz;
        let Vx0 = this.Vx0;
        let Vy0  = this.Vy0;
        let Vz0 = this.Vz0;
        let s = this.s;
        let density = this.density;

        diffuse(1, Vx0, Vx, visc, dt);
        diffuse(2, Vy0, Vy, visc, dt);
        diffuse(3, Vz0, Vz, visc, dt);

        project(Vx0, Vy0, Vx, Vy);

        advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt);
        advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt);
        advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt);

        project(Vx, Vy, Vz, Vx0, Vy0);

        diffuse(0, s, density, diff, dt);
        advect(0, density, s, Vx, Vy, Vz, dt);
    }

    // method to add density
    addDensity(x, y, z, amount) {
        let index = IX(x, y, z);
        this.density[index] += amount;
    }

    // method to add velocity
    addVelocity(x, y, z, amountX, amountY, amountZ) {
        let index = IX(x, y, z);
        this.Vx[index] += amountX;
        this.Vy[index] += amountY;
        this.Vz[index] += amountZ;
    }

    // function to render density
    renderD() {
        colorMode(HSB, 255);
    
        for (let i = 0; i < N; i++) {
            for (let j = 0; j < N; j++) {
                for (let k = 0; k < N; k++){
                    let x = i * SCALE;
                    let y = j * SCALE;
                    let z = k * SCALE;
                    let d = this.density[IX(i, j, k)];
                    fill((d + 50) % 255,200,d);
                    noStroke();
                    box(x, y, z, SCALE);
                }
            }
        }
    }
    
    // function to render velocity
    renderV() {

        for (let i = 0; i < N; i++) {
            for (let j = 0; j < N; j++) {
                for (let k = 0; k < N; k++) {
                    let x = i * SCALE;
                    let y = j * SCALE;
                    let z = k * SCALE;
                    let vx = this.Vx[IX(i, j, k)];
                    let vy = this.Vy[IX(i, j, k)];
                    let vz = this.Vz[IX(i, j, k)]
            
                    if (!(abs(vx) < 0.1 && abs(vy) <= 0.1)) {
                        beginShape(LINES);
                        stroke(255);
                        curveVertex(x, y, z);
                        curveVertex(x+vx*SCALE, y+vy*SCALE, z+vz*SCALE);
                        endShape();
                    }
                }
            }
        }
      }
}

