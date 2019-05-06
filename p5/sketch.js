/// <reference path="./p5.global-mode.d.ts" />

let fluid = new Fluid();

function settings() {
  size(N*SCALE, N*SCALE);
}

function setup() {
    createCanvas(1200, 1200, WEBGL);
    canvas = createGraphics(100, 100);
    canvas.background(0);
    translate(600, 350);
    fluid = new Fluid(0.2, 0, 0.0000001);
}

function draw() {
    background(255);
    let cx = int(0.5*width/SCALE);
    let cy = int(0.5*height/SCALE);
    for (let k = -1; k <= 1; k++){
        for (let i = -1; i <= 1; i++) {
            for (let j = -1; j <= 1; j++) {
                fluid.addDensity(cx+i, cy+j, k, random(50, 150));
        }
    }
    }
    
    for (let i = 0; i < 2; i++) {
        let angle = noise(t, t, t) * TWO_PI * 2;
        let v = p5.Vector.fromAngle(angle);
        v.mult(0.2);
        t += 0.01;
        fluid.addVelocity(cx, cy, 1, v.x, v.y, v.z);
    }

    fluid.step();
    fluid.renderD();
}