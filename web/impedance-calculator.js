function calcImpedance() {
    const m_to_ft = (m) => m * 3.28084;
    const ft_to_m = (ft) => ft * 0.30479999;
    const mi_to_m = (mi) => mi * 1.60934 * 1000;
    const km_to_mi = (km) => km * 0.621371;
    const dist = (x1, y1, x2, y2) => Math.sqrt(Math.pow(x1 - x2, 2) + Math.pow(y1 - y2, 2));

    // System parameters
    const f = parseFloat(document.getElementById("frequency").value); // Hz

    // Conductor positions (feet)
    //const apos = [parseFloat(document.getElementById("ax".value)), parseFloat(document.getElementById("ay").value)];
    const apos = [parseFloat(document.getElementById("ax").value), parseFloat(document.getElementById("ay").value)];
    const bpos = [parseFloat(document.getElementById("bx").value), parseFloat(document.getElementById("by").value)];
    const cpos = [parseFloat(document.getElementById("cx").value), parseFloat(document.getElementById("cy").value)];
    const n1pos = [parseFloat(document.getElementById("n1x").value), parseFloat(document.getElementById("n1y").value)];
    const n2pos = [parseFloat(document.getElementById("n2x").value), parseFloat(document.getElementById("n2y").value)];
    
    // Conductor properties
    const Ra = parseFloat(document.getElementById("Ra").value);           // Ohms per mile per conductor
    const Rb = Ra;
    const Rc = Ra;
    const GMRa = parseFloat(document.getElementById("GMRa").value);         // Geometric mean radius (feet)
    const GMRb = GMRa;
    const GMRc = GMRa;
    const ODa = parseFloat(document.getElementById("ODa").value);
    const bundle = parseFloat(document.getElementById("bundle").value);            // Conductors in bundle
    const bundle_spacing = parseFloat(document.getElementById("bundle_spacing").value);  // feet

    // Neutral conductor properties
    const Rn1 = parseFloat(document.getElementById("Rn1").value);            // Ohms per km
    const Rn2 = Rn1;
    const GMRn1 = parseFloat(document.getElementById("GMRn1").value);        // feet
    const ODn = parseFloat(document.getElementById("ODn").value);

    // values for dry earth
    Dkkp = 1
    if(document.getElementById("material").value == "seawater") {
        Dkkp = m_to_ft(38.25);  // feet
    } else if (document.getElementById("material").value == "swampyground") {
        Dkkp = m_to_ft(560);
    } else if (document.getElementById("material").value == "dampearth") {
        Dkkp = m_to_ft(850);
    } else if (document.getElementById("material").value == "dryearth") {
        Dkkp = m_to_ft(2690);
    } else if (document.getElementById("material").value == "pureslate") {
        Dkkp = m_to_ft(269000);
    } else if (document.getElementById("material").value == "sandstone") {
        Dkkp = m_to_ft(2690000);
    } else {
        Dkkp = m_to_ft(658.5*math.pow((parseFloat(document.getElementById("resistivity").value)/f),0.5))
    }

    const Rkp = 9.869E-7 * f;
    const Rk = (Ra / mi_to_m(1)) / bundle; // convert to ohms per meter per phase
    const Rn = (Rn1 / 1000);                  // convert to ohms per meter

    Dkkc = Math.pow(GMRa * Math.pow(bundle_spacing, (bundle - 1)), (1 / bundle));
    if (bundle >= 4){
        Dkkc = Dkkc * 1.091;
    }
    const Dkkg = GMRn1;

    const Zaa = Zbb = Zcc = math.add(Rk + Rkp, math.multiply(math.multiply(math.complex(0, 1),(2 * Math.PI * 60)*(2E-7)),(math.log(Dkkp / Dkkc))));
    const Zn1n1 = Zn2n2 = math.add( Rn + Rkp, math.multiply(math.complex(0, 1),(2 * Math.PI * 60)*(2E-7)*(math.log(Dkkp / GMRn1))));

    const Zab = Zbc = math.add(Rkp, math.multiply(math.complex(0, 1),(2 * Math.PI * 60)*(2E-7)*(math.log(dist(apos[0], apos[1], bpos[0], bpos[1] - Dkkp) / dist(apos[0], apos[1], bpos[0], bpos[1])))));
    const Zac = math.add(Rkp, math.multiply(math.complex(0, 1),(2 * Math.PI * 60)*(2E-7)*(math.log(dist(apos[0], apos[1], cpos[0], cpos[1] - Dkkp) / dist(apos[0], apos[1], cpos[0], cpos[1])))));

    const Zan1 = Zcn2 = math.add(Rkp, math.multiply(math.complex(0, 1),(2 * Math.PI * 60)*(2E-7)*(math.log(dist(apos[0], apos[1], n1pos[0], n1pos[1] - Dkkp) / dist(apos[0], apos[1], n1pos[0], n1pos[1])))));
    const Zan2 = Zcn1 = math.add(Rkp, math.multiply(math.complex(0, 1),(2 * Math.PI * 60)*(2E-7)*(math.log(dist(apos[0], apos[1], n2pos[0], n2pos[1] - Dkkp) / dist(apos[0], apos[1], n2pos[0], n2pos[1])))));
    const Zbn1 = Zbn2 = math.add(Rkp, math.multiply(math.complex(0, 1),(2 * Math.PI * 60)*(2E-7)*(math.log(dist(bpos[0], bpos[1], n1pos[0], n1pos[1] - Dkkp) / dist(bpos[0], bpos[1], n1pos[0], n1pos[1])))));

    const Zn1n2 = math.add(Rkp, math.multiply(math.complex(0, 1),(2 * Math.PI * 60)*(2E-7)*(math.log(dist(n1pos[0], n1pos[1], n2pos[0], n2pos[1] - Dkkp) / dist(n1pos[0], n1pos[1], n2pos[0], n2pos[1])))));
    console.log(Zan2);

    const Za = [
        [Zaa, Zab, Zac],
        [Zab, Zbb, Zbc],
        [Zac, Zbc, Zcc]
    ];

    const Zb = [
        [Zan1, Zan2],
        [Zbn1, Zbn2],
        [Zcn1, Zcn2]
    ];

    const Zc = [
        [Zan1, Zbn1, Zcn1],
        [Zan2, Zbn2, Zcn2]
    ];

    const Zd = [
        [Zn1n1, Zn1n2],
        [Zn1n2, Zn2n2]
    ];

    const Zp = math.subtract(Za, math.multiply(math.multiply(Zb, math.inv(Zd)), Zc));
    
    const Zs = math.divide(math.add( math.add(Zp[0][0], Zp[1][1]), Zp[2][2]),3);
    const Zm = math.divide(math.add( math.add(Zp[0][1], Zp[0][2]), Zp[1][2]),3);

    console.log(Zs);
    console.log(Zm);
    const Zpsym = [
        [Zs, Zm, Zm],
        [Zm, Zs, Zm],
        [Zm, Zm, Zs]
    ];

    const a = math.complex(math.cos(120 * (Math.PI / 180)), math.sin(120 * (Math.PI / 180)));
    const Ainv = [
        [1, 1, 1],
        [1, a, math.pow(a, 2)],
        [1, math.pow(a, 2), a]
    ];
    const A = [
        [1, 1, 1],
        [1, math.pow(a, 2), a],
        [1, a, math.pow(a, 2)]
    ];

    const Zshat = math.multiply((1 / 3), math.multiply(math.multiply(Ainv, Zpsym), A));  // ohms per meter
    console.log(`Z0 = ${Zshat[0][0].re * 1000} + j${Zshat[0][0].im * 1000} ohm/km (zero sequence)`);
    console.log(`Z1 = ${Zshat[1][1].re * 1000} + j${Zshat[1][1].im * 1000} ohm/km (positive sequence)`);
    console.log(`Z2 = ${Zshat[2][2].re * 1000} + j${Zshat[2][2].im * 1000} ohm/km (negative sequence)`);

    document.getElementById("Z0").textContent = `${Zshat[0][0].re * 1000} + j${Zshat[0][0].im * 1000} ohm/km (zero sequence)`;
    document.getElementById("Z1").textContent = `${Zshat[1][1].re * 1000} + j${Zshat[1][1].im * 1000} ohm/km (positive sequence)`;
    document.getElementById("Z2").textContent = `${Zshat[2][2].re * 1000} + j${Zshat[2][2].im * 1000} ohm/km (negative sequence)`;


    const ra = ODa / (2 * 12); // inches to feet
    const rn = ODn / (2 * 12);

    let Dsc = Math.pow(ra * Math.pow(bundle_spacing, (bundle - 1)), (1 / bundle));
    if (bundle >= 4) {
        Dkkc *= 1.091;
    }
    let Dscn = rn;

    const Paa = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((apos[1] * 2) / Dsc);
    const Pbb = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((bpos[1] * 2) / Dsc);
    const Pcc = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((cpos[1] * 2) / Dsc);
    const Pn1n1 = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((n1pos[1] * 2) / Dscn);
    const Pn2n2 = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((n2pos[1] * 2) / Dscn);

    const Pab = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(apos[0], apos[1], bpos[0], (-1 * bpos[1]))) / (dist(apos[0], apos[1], bpos[0], bpos[1])));
    const Pac = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(apos[0], apos[1], cpos[0], (-1 * cpos[1]))) / (dist(apos[0], apos[1], cpos[0], cpos[1])));
    const Pbc = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(bpos[0], bpos[1], cpos[0], (-1 * cpos[1]))) / (dist(bpos[0], bpos[1], cpos[0], cpos[1])));

    const Pan1 = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(apos[0], apos[1], n1pos[0], (-1 * n1pos[1]))) / (dist(apos[0], apos[1], n1pos[0], n1pos[1])));
    const Pan2 = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(apos[0], apos[1], n2pos[0], (-1 * n2pos[1]))) / (dist(apos[0], apos[1], n2pos[0], n2pos[1])));
    const Pbn1 = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(bpos[0], bpos[1], n1pos[0], (-1 * n1pos[1]))) / (dist(bpos[0], bpos[1], n1pos[0], n1pos[1])));
    const Pbn2 = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(bpos[0], bpos[1], n2pos[0], (-1 * n2pos[1]))) / (dist(bpos[0], bpos[1], n2pos[0], n2pos[1])));
    const Pcn1 = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(cpos[0], cpos[1], n1pos[0], (-1 * n1pos[1]))) / (dist(cpos[0], cpos[1], n1pos[0], n1pos[1])));
    const Pcn2 = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(cpos[0], cpos[1], n2pos[0], (-1 * n2pos[1]))) / (dist(cpos[0], cpos[1], n2pos[0], n2pos[1])));

    const Pn1n2 = (1 / (2 * Math.PI * 8.854E-12)) * Math.log((dist(n1pos[0], n1pos[1], n2pos[0], (-1 * n2pos[1]))) / (dist(n1pos[0], n1pos[1], n2pos[0], n2pos[1])));

    const Pa = [
        [Paa, Pab, Pac],
        [Pab, Pbb, Pbc],
        [Pac, Pbc, Pcc]
    ];

    const Pb = [
        [Pan1, Pan2],
        [Pbn1, Pbn2],
        [Pcn1, Pcn2]
    ];

    const Pc = [
        [Pan1, Pbn1, Pcn1],
        [Pan2, Pbn2, Pcn2]
    ];

    const Pd = [
        [Pn1n1, Pn1n2],
        [Pn1n2, Pn2n2]
    ];
        
    const Cp = math.inv(math.subtract(Pa, math.multiply(math.multiply(Pb, math.inv(Pd)), Pc)));

    const Caa = math.divide(math.add( math.add(Cp[0][0], Cp[1][1]), Cp[2][2]),3);
    const Cab = math.divide(math.add( math.add(Cp[0][1], Cp[0][2]), Cp[1][2]),3);

    const Cpsym = [
        [Caa, Cab, Cab],
        [Cab, Caa, Cab],
        [Cab, Cab, Caa]
    ];

    //Cphat = np.multiply(1j*(2*np.pi*f)*1000,(1/3)*np.matmul(np.matmul(Ainv,Cpsym),A))
    // Cshat = (1/3)*np.matmul(np.matmul(Ainv,Cpsym),A)  //ohms per meter
    // 1j*2*np.pi*Cpsym
    const Cphat = math.multiply(math.multiply(math.complex(0,1), (2*Math.PI*f/3)), math.multiply(math.multiply(Ainv, Cpsym), A)); //math.multiply(math.multiply(math.complex(0, 1),(2 * Math.PI * f * 1000)), math.multiply((1 / 3), math.multiply(math.multiply(Ainv, Cpsym), A)));
    console.log(Cphat);
    console.log(`Y0 = ${Cphat[0][0].im * 1000} ohm/km (zero sequence)`);
    console.log(`Y1 = ${Cphat[1][1].im * 1000} ohm/km (positive sequence)`);
    console.log(`Y2 = ${Cphat[2][2].im * 1000} ohm/km (negative sequence)`);

    document.getElementById("Y0").textContent = `j${Cphat[0][0].im * 1000} ohm/km (zero sequence)`;
    document.getElementById("Y1").textContent = `j${Cphat[1][1].im * 1000} ohm/km (positive sequence)`;
    document.getElementById("Y2").textContent = `j${Cphat[2][2].im * 1000} ohm/km (negative sequence)`;

}

function dropdownChangeHandler() {
    if(document.getElementById("material").value=="custom") {
        document.getElementById("resistivitylabel").style="display:inline;";
        document.getElementById("resistivity").style="display:inline;";
        document.getElementById("resistivityunits").style="display:inline;";
    } else {
        document.getElementById("resistivitylabel").style="display:none;";
        document.getElementById("resistivity").style="display:none;";
        document.getElementById("resistivityunits").style="display:none;";
    }
}
