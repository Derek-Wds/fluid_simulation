/// <reference path="./p5.global-mode.d.ts" />

function setup() {
    createCanvas(1200, 1200);
    canvas = createGraphics(100, 100);
    canvas.background(0)
    translate(600, 350);
}

function draw() {
    background(255);
    let jack = new Person('jack');
    textSize(32);
    text(jack.name, 10, 30);
}

class Person {
    constructor(name) {
        this.name = name;
    }
}