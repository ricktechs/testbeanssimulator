'''
================================================================================

MODELLING AND OPTIMIZATION OF DESIGN PARAMETERS AND\
PERFORMANCE OF A STATIONARY SOLAR THRESHER FOR COMMON BEAN
            (PHASEOLUS VULGARIS L)

            PYTHON 3.8 v2
================================================================================
'''

# libraries and modules
from flask import Flask
from flask import render_template, request
import datetime
import time
import math
import os
import os.path
from os import system
app = Flask(__name__)


# global constants
g = 10
kb = 0.4
u = 0.375
ks = 2.53
#kt = 9.14*pow(10,-4)
kt = 0.0021
km = 0.69
ke = 0.91
kn = 2.25
kf = 0.06
kd = 7.94*pow(10,-3)
Z = 0.7472
ka = 0.15
html_output = []
form_fields = {}
errors = {}
random_id = {}


@app.route('/model', methods=['GET', 'POST'])
def main():

    try:
        errors.clear()
        if request.method == 'GET':
            html_output.clear()
            # form_fields.clear()
            # print(time.mktime(datetime.datetime.now().timetuple()) * 1000)
        elif request.method == 'POST':
            vb = float(request.form['vb'])
            lc = float(request.form['lc'])
            fr = float(request.form['fr'])
            mc = float(request.form['mc'])
            bd = float(request.form['bd'])
            wc = float(request.form['wc'])
            a1 = float(request.form['a1'])
            a2 = float(request.form['a2'])
            b1 = float(request.form['b1'])
            b2 = float(request.form['b2'])
            d1 = float(request.form['d1'])
            D = float(request.form['D'])
            n = float(request.form['n'])
            y = float(request.form['y'])
            c = float(request.form['c'])
            b = float(request.form['b'])/100

            # update random id
            random_id.update(
                {"current": time.mktime(datetime.datetime.now().timetuple()) * 1000})

            # update form fields
            form_fields.update(
                {"vb": vb, "lc": lc, "fr": fr, "mc": mc, "bd": bd, "wc": wc, "a1": a1, "a2": a2, "b1": b1, "b2": b2, "d1": d1, "D": D, "n": n, "y": y, "c": c, "b": b})
            vg = (2/3)*vb
            vd = vg
            bw = bd/(1-(0.01*b))
            # Velocity of the crop in the threshing zone
            vc = kb*vb

            # Dwell time
            td = (1/kb)*(lc/vb)

            # Crop Stream Thickness
            S = ks*math.sqrt((fr/(vb*bw)))

            # Force Analysis of the crop stream in the threshing zone
            F = ka*fr*vb

            # b is the moisture content of the wet crop

            ##########Power requirement model###########
            # Power to detach the grains from the panicles
            #E = ke*(pow((fr*vb),3)*(1-b))/(bd*pow(lc,2))
            E = ke*(math.sqrt(vb)*pow(fr, 1.5))/math.sqrt(bw)
            p_i = E/td

            #Also, P_i = kb*(E*vc/lc)
            # Power to overcome frictional force
            pf = kf*fr*pow(vb, 2)

            # Power to turn unloaded cylinder
            pr = ((44/7)*mc*y*n*(g+(2*pow(vb,2)/D))/60000)

            # Total Power
            pt = p_i + pr + pf

            # Modelling threshing process
            # Threshing frequency
            l1 = kt*((pow(vb,2)*bd*D)/((1-b)*fr))

            #l1 = (pow(vb,2)*bd*D)/((1-b)*fr)

            # Velocity of grain after impact
            #fc = m*(vb-vg)

            # Power efficiency
            e = 2*(vb-vg)/vg

            # Threshing efficiency

            # Dwell time in the thresher (tdt)
            tdt = 3*lc/(2*vb)
            #te = 1-math.exp((-l1*td))
            ktl = kt*-1
            te = 1 - math.exp((ktl*bd*D*vb*lc)/((1-b)*fr))
            # Threshing loss
            tl = 1 - te

            tl2 = math.exp((-1.5*kt*bd*D*vb*lc)/((1-b)*fr))

            # Demage Model
            #ld = kd*math.sqrt((pow(D,2)*vb*bw/fr))*math.sqrt((pow(vd,3)*bw/fr))*vb/vd
            #gd = math.exp(-0.5*kd*(1-b)*wc*vb*lc/fr)
            gd = math.exp((-0.46*kd*bd*D*vb*lc)/((1-b)*fr))
            # demage fraction
            #df = math.exp((-ld*td))

            # total grain loss
            #tgl = math.exp((-0.5*kt*(bd*(1-b)*vb*wc*D*lc))/(c*fr))
            tgl = tl + gd

            # Grain migration parameter
            l2 = (1/kn)*math.sqrt((g+(2*pow(vb,2)/D)/math.sqrt((fr*(1-b)/(bd*vb)))))

            # probability of grain passage
            #P = (l3*3*b1)/(2*vb)
            P = (a1-a2-d1)*(b1-b2-d1)/(a1*b1)

            # number of grains passing through the concave opening in one second
            l3 = (2*vb*(a1-a2-d1)*(b1-b2-d1))/(3*a1*pow(b1,2))

            # separation efficiency
            Se = 1-(((l1*l3*(l3-l1)*math.exp((-l2*td)))+(l2*l1*(l1-l2)*math.exp((-l3*td)))+ (l2*l3*(l2-l3)*math.exp((-l1*td))))/((l1-l2)*(l3-l2)*(l3-l1)))

            # output capacity
            ct = km * fr * Z * Se*3600

            # Outputs
            output = {"pf": pf, "pr": pr, "pt": pt,"te": te,"tl": tl, "gd": gd, "tgl": tgl, "l3": l3,"ct": ct}
            output.update(form_fields)
            html_output.append(output)

    except:
        errors.update(
            {"message": "please fill all fields correctly and try again"})

    return render_template('model1.html', random_id=random_id, title='calculate', errors=errors, heading="Beans Simulation Model", result=html_output, form_fields=form_fields)
    main()


@app.route('/')
def modeler(name=None):
    return render_template('index.html', name=name)


@app.route('/about')
@app.route('/about/')
def about(name=None):
    return render_template('about.html', name=name)


@app.route('/services')
def services(name=None):
    return render_template('services.html', name=name)


@app.route('/contact')
def contact(name=None):
    return render_template('contact.html', name=name)


@app.route('/opt')
def opt(name=None):
    return render_template('opt.html', name=name)
