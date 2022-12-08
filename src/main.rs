use plotters::prelude::*;
extern crate nalgebra as na;
use fastrand;

struct Quadtree {
    center: na::Vector2<f64>,
    halfsize: f64,
    mass: f64,
    com: na::Vector2<f64>,
    children: Option<Vec<Quadtree>>,
    bodies: Option<Vec<(na::Vector2<f64>, f64)>>,
}
impl Quadtree {
    //const LIMIT: u32 = 1;
    fn new(center: na::Vector2<f64>, halfsize: f64) -> Quadtree {
        Quadtree {
            center: center,
            halfsize: halfsize,
            mass: 0.0,
            com: center,
            children: None,
            bodies: None,
        }
    }
    fn generate(points: &Vec<na::Vector2<f64>>, masses: &Vec<f64>) -> Quadtree {
        let mut minx = f64::MAX;
        let mut maxx = f64::MIN;
        let mut miny = f64::MAX;
        let mut maxy = f64::MIN;
        let mut com = na::vector![0.0, 0.0];
        for i in 0..points.len() {
            if points[i][0] < minx {
                minx = points[i][0]
            }
            if points[i][0] > maxx {
                maxx = points[i][0]
            }
            if points[i][1] < miny {
                miny = points[i][1]
            }
            if points[i][1] > maxy {
                maxy = points[i][1]
            }

            com += points[i]; // * masses[i];
        }
        //TODO: set abs limit for quad bounds

        com = com / points.len() as f64;
        let dif = (maxx - com[0])
            .max(maxy - com[1])
            .max(com[0] - minx)
            .max(com[1] - miny);
        let mut tree = Quadtree::new(com, dif * 1.1);
        for i in 0..points.len() {
            tree.insert(points[i], masses[i]);
        }
        tree
    }
    fn subdivide(&mut self) {
        let quartersize = self.halfsize / 2.0;
        let c1: na::Vector2<f64> = self.center + na::Vector2::new(quartersize, quartersize);
        let c2: na::Vector2<f64> = self.center + na::Vector2::new(-quartersize, quartersize);
        let c3: na::Vector2<f64> = self.center + na::Vector2::new(-quartersize, -quartersize);
        let c4: na::Vector2<f64> = self.center + na::Vector2::new(quartersize, -quartersize);

        self.children = Some(vec![
            Quadtree::new(c1, quartersize),
            Quadtree::new(c2, quartersize),
            Quadtree::new(c3, quartersize),
            Quadtree::new(c4, quartersize),
        ])
    }
    fn containspoint(&self, point: na::Vector2<f64>) -> bool {
        //bias towards lower boundary
        return (self.center[0] - self.halfsize) <= point[0]
            && point[0] < (self.center[0] + self.halfsize)
            && (self.center[1] - self.halfsize) <= point[1]
            && point[1] < (self.center[1] + self.halfsize);
    }
    fn insert(&mut self, point: na::Vector2<f64>, mass: f64) {
        //updates center of mass of current quad
        let locsum = self.com * self.mass + point * mass;
        self.mass += mass;
        self.com = locsum / self.mass;

        match self.bodies {
            None => self.bodies = Some(vec![(point, mass)]),
            _ => self.bodies.as_mut().unwrap().push((point, mass)),
        }

        //quick check to see if quad is empty expect for new addition, if it isnt, subdivide
        if self.mass != mass {
            //subdivides in children arent already defined
            match self.children {
                None => self.subdivide(),
                _ => (),
            }
            for p in self.bodies.as_ref().unwrap() {
                for child in self.children.as_mut().unwrap() {
                    if child.containspoint(p.0) {
                        child.insert(p.0, p.1);
                        break;
                    }
                }
            }
            self.bodies = None;
        }
    }
    //check barnes hut variables
    fn barnes_hut(&self, point: na::Vector2<f64>, theta: f64, g: f64) -> na::Vector2<f64> {
        if self.halfsize / modulus_squared(self.com - point).sqrt() < theta
            || self.bodies != None
            || self.children.is_none() && point != self.com
        {
            //println!(
            //    "Theta = {}/{}",
            //    self.halfsize,
            //    modulus_squared(self.com - point).sqrt()
            //);
            acceleration_function(point, self.com, self.mass, g)
        } else {
            let mut temp = na::vector!(0.0, 0.0);
            for child in self.children.as_ref().unwrap() {
                if child.mass != 0.0 {
                    temp += child.barnes_hut(point, theta, g)
                }
            }
            temp
        }
    }

    fn draw_prep(&self) -> Vec<[(f64, f64); 2]> {
        let halfvec = na::vector![self.halfsize, self.halfsize];
        let min = self.center - halfvec;
        let max = self.center + halfvec;

        let rectangle: [(f64, f64); 2] = [(min[0], min[1]), (max[0], max[1])];

        let mut rectangles = vec![rectangle];
        match self.children {
            None => return rectangles,
            Some(_) => {
                for child in self.children.as_ref().unwrap() {
                    rectangles.append(&mut child.draw_prep())
                }
            }
        }
        rectangles
    }

    //*area.draw(&Rectangle::new(
    //    [(min, min), (max, max)],
    //    Into::<ShapeStyle>::into(&GREEN),

    fn print(&self, indentation: usize) {
        println!(
            "{:indent$}COM:{}{} MASS: {}",
            "",
            self.com[0],
            self.com[1],
            self.mass,
            indent = indentation * 3
        );
        match self.children {
            None => (),
            _ => {
                for child in self.children.as_ref().unwrap() {
                    child.print(indentation + 1);
                }
            }
        }
    }
}

fn modulus_squared(deltap: na::Vector2<f64>) -> f64 {
    deltap[0].powf(2.0) + deltap[1].powf(2.0) + 0.001
}
fn acceleration_function(
    p1: na::Vector2<f64>,
    p2: na::Vector2<f64>,
    mass2: f64,
    g: f64,
) -> na::Vector2<f64> {
    let deltap = p2 - p1;
    let distance2 = modulus_squared(deltap);

    let mag = mass2 * g / distance2;
    let temp = deltap / modulus_squared(deltap).sqrt();
    temp * mag
}

fn sum_accelerations(
    bodies: &Vec<na::Vector2<f64>>,
    theta: f64,
    tree: Quadtree,
    g: f64,
) -> Vec<na::Vector2<f64>> {
    let mut acc: Vec<na::Vector2<f64>> = vec![];
    for i in 0..bodies.len() {
        let body = bodies[i];
        acc.push(tree.barnes_hut(body, theta, g));
    }
    acc
}

fn euler_step(
    //this might be fucked up idk
    bodies: &Vec<na::Vector2<f64>>,
    masses: &Vec<f64>,
    velocities: &Vec<na::Vector2<f64>>,
    tstep: f64,
    theta: f64,
    g: f64,
) -> [Vec<na::Vector2<f64>>; 2] {
    let mut dv = vec![];
    let mut da = vec![];
    let tree = Quadtree::generate(bodies, masses);
    //tree.print(0);
    //println!("{:?}", bodies);
    let accelerations = sum_accelerations(bodies, theta, tree, g);
    for i in 0..bodies.len() {
        dv.push(velocities[i] * tstep);
        da.push(accelerations[i] * tstep);
    }

    return [dv, da];
}

fn runge_kutta(
    bodies: &mut Vec<na::Vector2<f64>>,
    masses: &Vec<f64>,
    velocities: &mut Vec<na::Vector2<f64>>,
    tstep: f64,
    theta: f64,
    g: f64,
) {
    let k1 = euler_step(&bodies, &masses, &velocities, tstep, theta, g);
    let mut bodies1 = vec![];
    let mut vel1 = vec![];
    for i in 0..bodies.len() {
        bodies1.push(bodies[i] + k1[0][i] / 2.0);
        vel1.push(velocities[i] + k1[1][i] / 2.0);
    }
    let k2 = euler_step(&bodies1, &masses, &vel1, tstep, theta, g);
    let mut bodies2 = vec![];
    let mut vel2 = vec![];
    for i in 0..bodies.len() {
        bodies2.push(bodies[i] + k2[0][i] / 2.0);
        vel2.push(velocities[i] + k2[1][i] / 2.0);
    }
    let k3 = euler_step(&bodies2, &masses, &vel2, tstep, theta, g);
    let mut bodies3 = vec![];
    let mut vel3 = vec![];
    for i in 0..bodies.len() {
        bodies3.push(bodies[i] + k3[0][i]);
        vel3.push(velocities[i] + k3[1][i]);
    }
    let k4 = euler_step(&bodies3, &masses, &vel3, tstep, theta, g);

    for i in 0..bodies.len() {
        bodies[i] += (k1[0][i] + 2.0 * k2[0][i] + 2.0 * k3[0][i] + k4[0][i]) / 6.0;
        velocities[i] += (k1[1][i] + 2.0 * k2[1][i] + 2.0 * k3[1][i] + k4[1][i]) / 6.0;
    }
}

fn draw_contour(
    bodies: &Vec<na::Vector2<f64>>,
    masses: &Vec<f64>,
    g: f64,
    theta: f64,
    minx: f64,
    maxx: f64,
    miny: f64,
    maxy: f64,
    contour: bool,
) {
    let root = BitMapBackend::new("images/coutour.jpg", (1000, 1000)).into_drawing_area();

    root.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root)
        .margin(20)
        .x_label_area_size(10)
        .y_label_area_size(10)
        .build_cartesian_2d(minx..maxx, miny..maxy)
        .unwrap();

    chart
        .configure_mesh()
        .disable_x_mesh()
        .disable_y_mesh()
        .draw()
        .unwrap();

    let plotting_area = chart.plotting_area();

    let range = plotting_area.get_pixel_range();

    let (pw, ph) = (range.0.end - range.0.start, range.1.end - range.1.start);
    let (xr, yr) = (chart.x_range(), chart.y_range());

    let step = (
        (xr.end - xr.start) / pw as f64,
        ((yr.end - yr.start) / ph as f64),
    );

    let tree = Quadtree::generate(bodies, masses);

    for x in 0..pw {
        for y in 0..ph {
            let temp = modulus_squared(tree.barnes_hut(
                na::vector![xr.start + step.0 * x as f64, yr.start + step.1 * y as f64],
                theta,
                g,
            ))
            .sqrt();
            let color = {
                match contour {
                    false => HSLColor(0.66, 0.0, temp.sqrt() / 7.0),
                    true => HSLColor(0.66, 0.0, (temp.sqrt() * 6.0).round() / 42.0),
                }
            };
            plotting_area
                .draw_pixel(
                    (xr.start + step.0 * x as f64, yr.start + step.1 * y as f64),
                    &color, //temp.ln() * 5.0).round() / 50.0, (temp.sqrt() * 6.0).round() / (36.0)
                )
                .unwrap();
        }
    }
}

fn draw_bodies(
    bodyhist: Vec<Vec<na::Vector2<f64>>>,
    dt: f64,
    minx: f64,
    maxx: f64,
    miny: f64,
    maxy: f64,
    masses: Vec<f64>,
    drawbox: bool,
    res: i32,
) //dt is ms per frame
{
    //fuck, data type might be incompatible
    let area = BitMapBackend::gif(
        "images/animated.gif",
        (1000, 1000),
        (dt * 1000.0).round() as u32, /* Each frame show 1s */
    )
    .unwrap()
    .into_drawing_area();
    let length = bodyhist.len();
    for i in 0..length {
        //println!("{} out of {}", i, length);
        area.fill(&BLACK).unwrap();

        let mut ctx = ChartBuilder::on(&area)
            .build_cartesian_2d::<std::ops::Range<f64>, std::ops::Range<f64>>(
                minx..maxx,
                miny..maxy,
            )
            .unwrap();

        ctx.draw_series(
            bodyhist[i]
                .iter()
                .map(|point| Circle::new((point[0], point[1]), 1, WHITE.filled())),
        )
        .unwrap();

        if drawbox {
            let tree = Quadtree::generate(&bodyhist[i], &masses);
            let rects = tree.draw_prep();
            ctx.draw_series(
                rects
                    .iter()
                    .map(|rect: &[(f64, f64); 2]| Rectangle::new(*rect, &BLUE)),
            )
            .unwrap();
        }

        area.present().unwrap();
    }
}

fn cloud(
    points: &mut Vec<na::Vector2<f64>>,
    velocities: &mut Vec<na::Vector2<f64>>,
    masses: &mut Vec<f64>,
    center_index: usize,
    n: i32,
    r_avd: f64,
    r_div: f64,
    v_var: f64,
    max_mass_ratio: f64,
    g: f64,
) {
    let rng = fastrand::Rng::new();

    for _i in 0..n {
        let r = r_avd + (rng.f64() - 0.5) * 2.0 * r_div;
        let theta = rng.f64() * 2.0 * std::f64::consts::PI;
        let v_temp = ((masses[center_index] + n as f64 * max_mass_ratio / 0.5) * g / r).sqrt();
        let v_base = v_temp + v_temp * v_var * (rng.f64() - 0.5);

        points.push(na::vector![
            points[center_index][0] + r * theta.cos(),
            points[center_index][1] + r * theta.sin()
        ]);
        velocities.push(na::vector![
            velocities[center_index][0] + v_base * (theta + std::f64::consts::PI / 2.0).cos(),
            velocities[center_index][1] + v_base * (theta + std::f64::consts::PI / 2.0).sin()
        ]);
        masses.push(rng.f64() * max_mass_ratio * masses[center_index]);
    }
}
fn main() {
    let tstep = 0.01; //ammount of time per simulation step
    let t: f64 = 10.0; //total simulation time
    let theta = 0.25; //ratio of half-sidelength / distance used in BH alg, lower values are more accurate, 0 is direct force summation
    let g = 0.1; //gravitational constant
    let playback = 5; //changes speed of simulation
    let halfsize = 10.0; //size in each direction displayed
    let drawbox = true; //decides whether or not to draw quadtree in animation
    let contour = true; //decides whether or not to draw a contour plot or just direct pixle values
    let res = 750;

    let frames = (t / tstep).round() as i32;

    //initializes coordinates of all bodies in system
    let mut points = vec![
        na::Vector2::new(0.0, 0.0),
        na::vector![9.0, 0.0],
        na::vector![-3.0, 0.0],
    ];
    //masses of all bodies
    let mut masses = vec![30.0, 20.0, 10.0, 0.01];
    // initial velocities
    let mut velocities = vec![
        na::vector![0.0, 0.0],
        na::vector![0.0, 0.5],
        na::vector![0.0, -0.8],
        na::vector![0.0, 0.0],
    ];
    //generate a bunch of light n bodies around object 1
    let n = 0;
    let r_avd = 6.0;
    let r_div = 0.3;
    let max_mass_ratio = 0.0002; // one part per million is similar to galaxy, ie five zeros
    cloud(
        &mut points,
        &mut velocities,
        &mut masses,
        0,
        n,
        r_avd,
        r_div,
        0.1,
        max_mass_ratio,
        g,
    );

    //funny normalization code, should put center of mass at origin and put average velocity = 0
    //interestingly, makes COM constant at (0,0) *assuming perfect simulation, which it isnt
    let mut com = na::vector![0.0, 0.0];
    let mut mass = 0.0;
    let mut cov = na::vector![0.0, 0.0];

    for i in 0..points.len() {
        com += points[i] * masses[i];
        cov += velocities[i] * masses[i];
        mass += masses[i];
    }
    cov /= mass;
    com /= mass;

    for i in 0..points.len() {
        points[i] -= com;
        velocities[i] -= cov;
    }
    points.push(na::vector![0.0, 0.0]);
    let momentum1 = modulus_squared(cov).sqrt() * mass;
    println!(
        "{} Bodies at {} Steps, Theta = {}",
        points.len(),
        frames,
        theta
    );
    let tic = std::time::Instant::now();
    //actual simulation code, records positions of bodies into a single vector
    let mut bodyhist = vec![points.to_vec()];
    for i in 0..frames {
        runge_kutta(&mut points, &masses, &mut velocities, tstep, theta, g);
        if i % playback == 0 {
            bodyhist.push(points.to_vec());
        }
    }

    println!("Runged in {} Seconds", tic.elapsed().as_secs_f64());

    let mut com = na::vector![0.0, 0.0];
    let mut mass = 0.0;
    let mut cov = na::vector![0.0, 0.0];
    for i in 0..points.len() {
        com += points[i] * masses[i];
        cov += velocities[i] * masses[i];
        mass += masses[i];
    }
    cov /= mass;
    com /= mass;
    let momentum2 = modulus_squared(cov).sqrt() * mass;
    println!(
        "COM Delta: {}   System Velocity Delta: {}   Average Momentum Change: {} Per Second",
        com,
        cov,
        (momentum2 - momentum1).abs() / t
    );

    let toc = std::time::Instant::now();
    draw_contour(
        &bodyhist[0],
        &masses,
        g,
        theta,
        -halfsize,
        halfsize,
        -halfsize,
        halfsize,
        contour,
    );

    draw_bodies(
        bodyhist,
        tstep * playback as f64,
        -halfsize,
        halfsize,
        -halfsize,
        halfsize,
        masses,
        drawbox,
        res,
    );

    println!("Rendered in {} Seconds", toc.elapsed().as_secs_f64());
}
