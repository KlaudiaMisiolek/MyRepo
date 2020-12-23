#!/usr/bin/env python3
import astroplan
from astroplan import Observer
import astropy
from astropy import coordinates, time
from astropy import units as u
import numpy as np
import matplotlib
from matplotlib import pylab as plt
from matplotlib.dates import DateFormatter, HourLocator, AutoDateLocator
from datetime import datetime, timedelta
import pandas as pd
from multiprocessing import Pool, cpu_count

def getCelestial(day, observer,observerLocation, target):
    searchWhich = 'nearest'
    target_rise = observer.target_rise_time(day, target, which=searchWhich)
    target_noon = observer.target_meridian_transit_time(day, target, which=searchWhich)
    target_set = observer.target_set_time(day, target, which=searchWhich)
    target_rise_deg = target.transform_to(astropy.coordinates.AltAz(obstime=target_rise, location=observerLocation)).az.degree
    target_noon_deg = target.transform_to(astropy.coordinates.AltAz(obstime=target_noon, location=observerLocation)).az.degree
    target_set_deg = target.transform_to(astropy.coordinates.AltAz(obstime=target_set, location=observerLocation)).az.degree
    # target_rise_diff = target_rise.datetime - day.datetime
    # print(matplotlib.dates.date2num(target_rise.datetime))
    target_rise = target_rise.datetime.replace(year=1970,month=1,day=1)
    target_noon = target_noon.datetime.replace(year=1970,month=1,day=1)
    target_set = target_set.datetime.replace(year=1970,month=1,day=1)
    return(
        {
            "day": day.datetime,
            "rise": target_rise, "noon": target_noon, "set": target_set,
            "rise_deg": target_rise_deg, "noon_deg": target_noon_deg, "set_deg": target_set_deg
        })
def getCelestialArgs(args):
    return(getCelestial(args["day"], args["observer"],args["observerLocation"], args["target"]))

def plotTarget(targetName, data):
    fig, ax = plt.subplots()
    locator = AutoDateLocator(minticks=12, maxticks=12)
    data.plot(x='day',y='rise',ax=ax)
    data.plot(x='day',y='set',ax=ax)
    data.plot(x='day',y='noon',ax=ax)
    ax.xaxis.set_major_formatter(DateFormatter('%m.%d'))
    ax.yaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.xaxis.set_major_locator(locator)
    plt.grid()
    plt.title(targetName + " time")
    plt.tight_layout()

    fig, ax = plt.subplots()
    data.plot(x='day',y='rise_deg',ax=ax)
    data.plot(x='day',y='set_deg',ax=ax)
    data.plot(x='day',y='noon_deg',ax=ax)
    ax.xaxis.set_major_formatter(DateFormatter('%m.%d'))
    # ax.yaxis.set_major_formatter(DateFormatter('%H:%M'))
    ax.xaxis.set_major_locator(locator)
    plt.grid()
    plt.title(targetName + " azimute")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    latitude = 50.06
    longitude = 19.94
    elevation = 220

    # daysOfTheYear = [datetime(2020, x,) for x in range(1,13)]
    # startDate = datetime(1678, 1,1)
    startDate = datetime(2020, 1,1,12)
    # endDate = datetime(2020, 1,20)
    # startDate = datetime(2020, 12,10)
    endDate = datetime(2020, 12,31)
    # endDate = datetime(2020, 12,31)
    __cores = cpu_count()
    daysOfTheYear = astropy.time.Time(pd.date_range(startDate,endDate,freq='d'))
    daysOfTheYear = daysOfTheYear[::6]

    observerLocation = astropy.coordinates.EarthLocation(lat=latitude*u.deg, lon=longitude*u.deg, height=elevation*u.m)
    # observer = astroplan.Observer(observerLocation,name="Krakow", timezone="Europe/Warsaw")
    observer = astroplan.Observer(observerLocation)
    sun = astropy.coordinates.get_sun(daysOfTheYear[0])
    moon =  astropy.coordinates.get_moon(daysOfTheYear[0])

    toPlot_sun = [{"day": day, "observer": observer, "observerLocation": observerLocation,
            "target": sun} for day in daysOfTheYear]
    toPlot_moon = [{"day": day, "observer": observer, "observerLocation": observerLocation,
            "target": moon} for day in daysOfTheYear]
    with Pool(processes=__cores) as pool:
        toPlot_sun = pool.map(getCelestialArgs, toPlot_sun)
        toPlot_moon = pool.map(getCelestialArgs, toPlot_moon)

    toPlot_sun = pd.DataFrame(toPlot_sun)
    toPlot_moon = pd.DataFrame(toPlot_moon)

    plotTarget("Sun",toPlot_sun)
    plotTarget("Moon",toPlot_moon)






