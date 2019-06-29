# (c) 2019 Anirban Banerjee
#Licensed under:
#GNU GENERAL PUBLIC LICENSE
#Version 3, 29 June 2007
from mpapbf import *

class astro ():
	pi = mpap('3.1415926535897932')
	pix2 = mpap('6.2831853071795865')
	D2R = pi/180
	R2D = mpap(180)/pi

	Longitude = mpap(81.889)
	Latitude = mpap(25.426)
	Elevation = mpap(0)

	Ls = mpap(0)
	Lm = mpap(0)
	Ms = mpap(0)
	Mm = mpap(0)

	month = ["January","February","March","April","May","June",
	   "July","August","September","October","November","December"]

	rashi = ["Mesha","Vrisha","Mithuna","Karka","Simha","Kanya","Tula",
	   "Vrischika","Dhanu","Makara","Kumbha","Meena"]

	day = ["Ravi","Soma","MAngal","Budh","Brihaspati","Shukra","Shani"]

	tithi = ["Prathamaa","Dvitiya","Tritiya","Chaturthi","Panchami",
			"Shashthi","Saptami","Ashtami","Navami","Dashami","Ekadashi",
			"Dvadashi","Trayodashi","Chaturdashi","Purnima","Pratipada",
			"Dvitiya","Tritiya","Chaturthi","Panchami","Shashthi",
			"Saptami","Ashtami","Navami","Dashami","Ekadashi","Dvadashi",
			"Trayodashi","Chaturdashi","Amaavasya"]

	karan = ["Bava","Baalava","Kaulava","Taitula","Garija","Vanija",
	   "Vishti","Shakuni","Chatushpada","Naga","Kimstughna"]

	yoga = ["Vishakumbha","Priti","Ayushman","Saubhagya","Shobhana",
	   "Atiganda","Sukarman","Dhriti","Shula","Ganda","Vriddhi",
	   "Dhruva","Vyaghata","Harshana","Vajra","Siddhi","Vyatipata",
	   "Variyan","Parigha","Shiva","Siddha","Saadhya","Shubha","Shukla",
	   "Brahma","Indra","Vaidhriti"]

	nakshatra = ["Ashvini","Bharani","Krittika","Rohini","Mrigashira","Ardra",
			"Punarvasu","Pushya","Ashlesa","Magha","Purva Phalguni","Uttara Phalguni",
			"Hasta","Chitra","Svaati","Vishakha","Anuradha","Jyeshtha","Mula",
			"Purva Ashadha","Uttara Ashadha","Shravana","Dhanishtha","Shatabhisha",
			"Purva Bhaadra","Uttara Bhaadra","Revati"]

	def calculate_ayanansha (self, d):
		t = (mpap(d) + 36523.5) / 36525
		o = mpap('259.183275') - mpap('1934.142008333206') * t + mpap('0.0020777778') * t * t
		L = mpap('279.696678') + mpap('36000.76892') *t + mpap('0.0003025') * t * t
		ayan = mpap(17.23) * (o * self.D2R).sin() + (L * self.D2R * 2).sin() * 1.27 - (mpap(5025.64) + mpap(1.11) * t) * t
		# Based on Lahiri
		ayan = (ayan-80861.27) / 3600
		return ayan

	def REV (self, x):
		return x - (x / 360).floor() * 360.0

	def lsun (self, d):
		#print ("d is ", d)
		w = mpap('282.9404') + mpap ('4.70935e-5') * d
		#print ("w is ", w)
		a = mpap(1.0)
		e = mpap(0.016709) - mpap('1.151e-9') * d
		M = self.REV(mpap('356.0470') + mpap('0.9856002585') * d)
		self.Ms = M
		self.Ls = w + M
		
		tmp = M * self.D2R
		E =  M + self.R2D * e * tmp.sin() * (e * tmp.cos() + 1)
		#print ("E is ", E)
		
		tmp = E * self.D2R
		x = tmp.cos() - e
		#print ("x is ", x)
		y = tmp.sin() * (mpap(1) - e * e).sqrt()
		#print ("y is ", y)
		
		r = (x * x + y * y).sqrt()
		#print ("r is ", r)
		#print ("y.atan2(x) is ", y.atan2(x))
		#print ("self.R2D * y.atan2(x) is ", self.R2D * y.atan2(x))
		v = self.REV (self.R2D * y.atan2(x))
		#print ("v is ", v)
		
		return self.REV(v+w)

	def lmoon (self, d):
		N = mpap('125.1228') - mpap('0.0529538083') * d
		i = mpap('5.1454')
		w = self.REV(mpap('318.0634') + mpap('0.1643573223') * d)
		a = mpap('60.2666')
		e = mpap('0.054900')
		M = self.REV(mpap('115.3654') + mpap('13.0649929509') * d)

		self.Mm = M
		self.Lm = N + w + M
		
		#Calculate Eccentricity anomaly
		tmp = M * self.D2R
		E = M + self.R2D * e * tmp.sin() * (e * tmp.cos() + 1)
		tmp = E * self.D2R
		Et = E - (E - self.R2D * e * tmp.sin() - M) / (mpap(1) - e * tmp.cos())
		
		while(E - Et > mpap('1e-3')):
			E = Et
			tmp = E * self.D2R
			Et = E - (E - self.R2D * e * tmp.sin() - M) / (mpap(1) - e * tmp.cos())
		
		tmp = E * self.D2R
		x = a * (tmp.cos() - e)
		y = a * (mpap(1) - e * e).sqrt() * tmp.sin()
		
		r = (x * x + y * y).sqrt()
		v = self.REV(self.R2D * y.atan2(x))
		
		tmp = self.D2R * N
		tmp1 = self.D2R * (v + w)
		tmp2 = self.D2R * i
		xec = r * (tmp.cos() * tmp1.cos() - tmp.sin() * tmp1.sin() * tmp2.cos())
		yec = r * (tmp.sin() * tmp1.cos() + tmp.cos() * tmp1.sin() * tmp2.cos())
		zec = r * tmp1.sin() * tmp2.sin()
		
		# Do some corrections
		D = self.Lm - self.Ls
		F = self.Lm - N
		
		lon = self.R2D * yec.atan2(xec) \
			- mpap('1.274')  * ((self.Mm - D * 2) * self.D2R).sin() \
			+ mpap('0.658')  * ((D * 2)*self.D2R).sin() \
			- mpap('0.186')  * (self.Ms * self.D2R).sin() \
			- mpap('0.059')  * ((self.Mm * 2 - D * 2) * self.D2R).sin() \
			- mpap('0.057')  * ((self.Mm - D * 2 + self.Ms) * self.D2R).sin() \
			+ mpap('0.053')  * ((self.Mm + D * 2) * self.D2R).sin() \
			+ mpap('0.046')  * ((D * 2 - self.Ms) * self.D2R).sin() \
			+ mpap('0.041')  * ((self.Mm - self.Ms) * self.D2R).sin() \
			- mpap('0.035')  * (D * self.D2R).sin() \
			- mpap('0.031')  * ((self.Mm + self.Ms) * self.D2R).sin() \
			- mpap('0.015')  * ((F * 2 - D * 2)*self.D2R).sin() \
			+ mpap('0.011')  * ((self.Mm - D * 4) * self.D2R).sin()
		return self.REV(lon)

	def calculate_panchanga (self, dd, mm, yy, hr, zhr):

		#Calculate day number since 2000 Jan 0.0 TDT
		d = mpap(367) * yy - mpap(7) * (mpap(yy) + mpap(mm + 9) / 12) / 4 + mpap(275) * mm / 9 + dd - 730530
		
		#Calculate Ayanamsa, moon and sun longitude
		ayanansha = self.calculate_ayanansha(d)
		#print ("ayanansha", ayanansha)
		vaara = self.day[int((d + 6 ) % 7)] #1 Jan 2000 was Saturday, so add 6; 
		d = d + (hr - zhr) / 24
		slon = self.lsun(d)
		mlon = self.lmoon (d)
		
		#print ("slon is ", slon)
		#print ("mlon is ", mlon)

		#Calculate Tithi and Paksha
		tmlon = mlon + (360 if (mlon < slon) else 0)
		tslon = slon
		n = int((tmlon - tslon) / 12)
		tithi = self.tithi[n]
		paksha = "Shukla" if (n <= 14) else "Krishna"
			
		#Calculate Nakshatra
		tmlon = self.REV(mlon + ayanansha)
		nakshatra = self.nakshatra[int(tmlon * 6 / 80)]
		
		#Calculate Yoga
		tmlon = mlon + ayanansha
		tslon = slon + ayanansha;
		yoga = self.yoga[int(self.REV(tmlon + tslon) * 6 / 80)]
		
		#Calculate Karana
		tmlon = mlon + 360 if (mlon < slon) else mpap(0)
		tslon = slon
		n = int((tmlon - tslon) / 6)
		if n <= 0:
			n = 10
		if n >= 57:
			n -= 50
		if n > 0 and n < 57:
			n = (n - 1) - int((n - 1) / 7) * 7
		karan = self.karan[n]
			
		#Calculate the rashi in which the moon is present
		tmlon = self.REV(mlon + ayanansha)
		rashi = self.rashi[int(tmlon / 30)]

		print(' Tithi: ' + tithi, '\n', paksha + ' paksha\n', vaara + 'vaara\n', 'Nakshatra: ' + nakshatra + '\n',\
				'Yoga: ' + yoga + '\n',\
				'Karana: ' + karan + '\n', 'Raashi: ' + rashi)

	def sunrise_equation (self, dd, mm, yy):
		#---day number since 2000 Jan 0.0 TDT
		d = mpap(367) * mpap(yy) - (mpap(7)) * (mpap(yy) + mpap(5001) + \
				(mpap(mm) - 9) // 7) // 4 + (mpap(275) * mm) // 9 + dd + 1729777
		#print ("Julian day to d is ", d)

		#d = Jdate is the Julian date
		#n is the number of days since Jan 1st, 2000 12:00.
		n = d - mpap('2451545.0') + mpap(68.184) / 86400

		#---mean solar noon
		#jSTAR is an approximation of mean solar time at d,
		#expressed as a Julian day with the day fraction.
		#self.Longitude is the longitude west (west is negative, 
		#east is positive) of the observer on the Earth in degrees
		jSTAR = n - self.Longitude/360.0

		#---solar mean anomaly
		M = ((mpap('357.5291') + mpap('0.98560028') * jSTAR) % 360) * self.D2R

		#---equation of the centre
		#C is the Equation of the center value needed to calculate lambda (see next equation).
		#1.9148 is the coefficient of the Equation of the Center for the planet the observer is on (in this case, Earth)
		C = (mpap('1.9148') * M.sin() + mpap('0.02') * (M * 2).sin() + mpap('0.0003') * (M*3).sin()) * self.D2R

		#---ecliptic longitude -- lambda
		#lambdaLong is the ecliptic longitude.
		#102.9372 is a value for the argument of perihelion.
		lambdaLong = (M + C + self.pi + mpap('102.9372') * self.D2R) % self.pix2

		#---solar transit
		#jTRANSIT is the Julian date for the local true solar transit (or solar noon).
		#2451545.0 is noon of the equivalent Julian year reference.
		#mpap(0.0053) * M.sin() - mpap('0.0069') * (lambdaLong * 2).sin() is a 
		#simplified version of the equation of time. The coefficients are fractional day minutes.
		jTRANSIT = mpap('2451545.0') + jSTAR + mpap(0.0053) * M.sin() - mpap('0.0069') * (lambdaLong * 2).sin()

		#---declination of the sun
		#delta is the declination of the sun
		#23.44° is Earth's maximum axial tilt toward the sun
		sineDelta = lambdaLong.sin() * (mpap('23.44') * self.D2R).sin()
		delta = sineDelta.asin()

		#---elevation correction (elevation is in metres)
		#This corrects for both apparent dip and terrestrial refraction. 
		#For example, for an observer at 10,000 feet, add (−115°/60°) or about −1.92° to −0.83°.
		elevationCorr = (mpap('-2.076') / 60) * self.Elevation.sqrt() * self.D2R
		
		#---hour angle
		#omega is the hour angle from the observer's zenith;
		#phi is the north latitude of the observer (north is positive, 
		#south is negative) on the Earth.
		phi = self.Latitude * self.D2R
		cosOmega = ((mpap('-0.83') * self.D2R + elevationCorr).sin() - phi.sin() * sineDelta) / \
					(phi.cos() * delta.cos())
		omega = cosOmega.acos()

		#calculate sunrise and sunset times
		jRISE = jTRANSIT - omega / self.pix2
		jSET  = jTRANSIT + omega / self.pix2

		#print ("Rise Time: ", jRISE)
		#print ("Set Time: ", jSET)

		jNOON = (jSET + jRISE) / 2
		dayLightHours = (jSET - jRISE) * 24
		#print ("jNOON is :", jNOON)
		riseToNoonHrs = (jNOON - jRISE) * 24

		noonToSetHrs = (jNOON - jRISE) * 24

		sriseHr = str( mpap(12) - riseToNoonHrs.ceil()) + ':' + ((mpap(1) - riseToNoonHrs.frac()) * 60).roundstr(2)
		ssetHr = str(int(noonToSetHrs) + 12) + ':' + (noonToSetHrs.frac() * 60).roundstr(2)
		print (" Sunrise at: ", sriseHr)
		print (" Sunset at: ", ssetHr)
		print (" Hours of daylight:", dayLightHours.roundstr(2))

		return
