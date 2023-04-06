import java.awt.Color;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import java.util.Random;

//**********************************************************************************************
// Author: Non-Euclidean Dreamer
// A physical Particle with all its basic properties
//**********************************************************************************************
public class Particle implements Comparable 
{
	static double g=6*Math.pow(10,-11),kc=9*Math.pow(10, 9),c=1,r;
	double[]loc, v;//for relativistic universes v is used for p
	double mass, radius,u;
	static double rulescope=Cluster.rulescope, metric=2;
	int  color;
	static boolean spherical;
	static int dim=2,background=Color.black.getRGB();
	String name;
	Random rand=new Random();
	double charge;
	public Particle(double[]l,double[]vel,double m, double r,int c)
	{
		loc=l.clone();
		v=vel.clone();
		mass=m;
		radius=r;
		dim=l.length;
		color=c;
		name="";
	}
	public Particle(double[] l, double[] vel, double m, double r, int[] is) 
	{
		loc=l.clone();
		v=vel.clone();
		mass=m;
		radius=r;
		dim=l.length;
		color=new Color(is[0],is[1],is[2]).getRGB();
		name="";
	}
	public Particle copy()
	{
		Particle obj= new Particle(loc.clone(),v.clone(),mass, radius,color);
		obj.name=name.concat("");
		return obj;
	}
	public static Particle zro()
	{
		return new Particle(new double[dim],new double[dim],10,10,background);
	}
	public static Particle proton(double[]loc, double[]v)
	{
		Particle out=new Particle(loc,v,1,1,Color.red.getRGB());
		out.setCharge(1);
		return out;
	}
	public void setV(double[]vel)
	{
		v=vel.clone();
	}
	public double distance(Particle obj,double norm)
	{
		double out=0;
		for(int i=0;i<dim;i++)
		{
			out+=Math.pow(Math.abs(loc[i]-obj.loc[i]), norm);
		}
		out=Math.pow(out, 1/norm);
		
		return out;
	}
	public double vdifference(Particle obj, double norm)
	{
		double out=0;
		for(int i=0;i<dim;i++)
		{
			out+=Math.pow(Math.abs(v[i]-obj.v[i]), norm);
		}
		out=Math.pow(out, 1/norm);
		
		return out;
	}
	
	//Changes the inertia&location
	public void pull(double[] acc,double t) 
	{
		System.out.println("acc="+acc[0]+","+acc[1]);
		for(int i=0;i<dim;i++)
		{
			v[i]+=t*acc[i];
			loc[i]+=t*v[i];
		}	
	}
	
	//To be used to merge two objects into one(when too close...)
	public boolean merge(Particle object,Cluster cluster, boolean zro) 
	{	//System.out.print("inertia:"+inertia[0]*mass+"+"+object.inertia[0]*object.mass);
		boolean out=false;
		double totalmass=mass+object.mass,absmass=Math.abs(mass)+Math.abs(object.mass);
		double u=0,ekin=0;
		if(!Cluster.fractal) {u=u()+object.u();
				ekin=ekin()+object.ekin();}
		else radius=Math.sqrt(mass*radius*radius+object.mass*object.radius*object.radius)/totalmass;
		charge+=object.charge;
		if(zro&&totalmass==0)
		{
			for(int i=0;i<dim;i++)
			{
			loc[i]=loc[i]+object.loc[i];
			loc[i]/=2;
			v[i]=v[i]+object.v[i];
			}
			color=Cluster.blur(color,object.color,0.5);
			out=true;
		}
		else if(zro)
		{
			for(int i=0;i<dim;i++)
			{
				loc[i]=object.loc[i];
				v[i]=v[i]+object.v[i]*(object.mass);
				v[i]/=totalmass;
			}
			
			color=Cluster.blur(color,object.color,Math.abs(mass/totalmass));
			
			mass=totalmass;
		}
		else if(totalmass!=0)
		{	
		for(int i=0;i<dim;i++)
		{
			if(spherical&&i==0&&(Math.abs(loc[i]-object.loc[i])>Math.PI)){if(loc[i]>0)loc[i]-=2*Math.PI;else loc[i]+=2*Math.PI;}
			loc[i]=loc[i]*Math.abs(mass)+object.loc[i]*Math.abs(object.mass);
			loc[i]/=absmass;
			if(spherical&&i==0) {loc[i]=(loc[i]+3*Math.PI)%(2*Math.PI)-Math.PI;}
			if(Cluster.relativistic)v[i]+=object.v[i];
			else 
			{
			v[i]=v[i]*(mass)+object.v[i]*(object.mass);
			v[i]/=totalmass;
			}
		}
		
		color=Cluster.blur(color,object.color,Math.abs(mass)/(Math.abs(mass)+Math.abs(object.mass)));
		
		mass=totalmass;
		}
		else if(mass==0)
		{
			for(int i=0;i<dim;i++)
			{
				loc[i]=loc[i]+object.loc[i];
				loc[i]/=2;
				v[i]=v[i]+object.v[i];
				v[i]/=2;
			}
			
			color=Cluster.blur(color,object.color,0.5);
			
		}
		else 
		{
			for(int i=0;i<dim;i++)
			{
				loc[i]=loc[i]*Math.abs(mass)+object.loc[i]*Math.abs(object.mass);
				loc[i]/=(Math.abs(mass)+Math.abs(object.mass));
				v[i]=v[i]*(mass)+object.v[i]*(object.mass);
			    
			}
			
			color=Cluster.blur(color,object.color,0.5);
			
			mass=totalmass;
			out=true;
		}//System.out.print("merge "+name+" with "+object.name);
	if(!Cluster.fractal)
	{	double vnorm=norm(v) ;
		u+=ekin-mass*c*(Math.sqrt(c*c+vnorm*vnorm)-c);
		 radius=Cluster.r(u,mass,Math.min(Math.PI,radius+object.radius));}
	
		
		if(name.equals("")&&object.name.equals("")) {if(Math.abs(mass)>Universe.ten) {System.out.print("new name=");name=randomSyllable();System.out.println(name);cluster.screens.add(new Screen( new BufferedImage(Screen.width,Screen.height,BufferedImage.TYPE_3BYTE_BGR),this));}else System.out.println("mass="+mass);}
		else 
		{
			if(name.equals("")&&!object.name.equals(""))cluster.replaceScreen(object,this);
			
			//we don't want the names to get too long://ToDo
			if(name.length()+object.name.length()>rand.nextDouble()*30)
			{
				//while(name.substring(name.length()-1).)
			}
			
			name=name.concat(object.name);
			//System.out.println("="+name);
		}//System.out.println("="+inertia[0]*mass);
		
		return out;
	}
	private double ekin() 
	{
		double vnorm=norm(v);
		return (Math.sqrt(c*c+vnorm*vnorm)-c)*mass*c;
	}
	private String randomSyllable() 
	{	
		final int cons=21;
		final String[] letters= {"p","b","m","f","v","t","d","n","th","dh","k","g","ng","ch","h","s","z","sh","zh","l","r","a","e","i","o","u","ae","oe","y","ai","au","ei","ja","je","jo","ju"};
		final int[][]standard= {{0,1,2,3,4},{5,6,7,8,9},{10,11,12,13,14}},
							hiss= {{15,16},{17,18}};
		final int[]special= {19,20},
							vowel= {21,22,23,24,25,26,27,28,29,30,31,32,33,34,35};
		final double[] startl= {0.25,0.4,0.3,0.05},//probability of i consonants at beginning
				twoconsstart={0.1 /*"pf"*/,0.1 /*"pr"*/, 0.2/*"ps"*/,
			 0.15/*"fr"*/,
			0.2 /*"sp"*/, 0.15  /*"sm"*/,0.1 /*"sr"*/},
				threeconsstart= {0.5/*"spr"*/,0.5 /*"spl"*/},
			 endl= {};
		int[]contains;
		final ArrayList<Integer> containsSt=new ArrayList<Integer>();
		containsSt.add(0);/*,con1,2,3,4,5,6,7,8,9,11,12);if (containsSt.*/
		String out="";
		
		//beginning:
		double start=rand.nextDouble();
		int clr = 0, voice,ssh = 0,lr = 0;
		if(start>=startl[0])
		{
			start-=startl[0];
			if(start>=startl[1])
			{
				start-=startl[1];
				if(start>=startl[2])
				{
					start-=startl[2];
					//3consonants
					{
						start/=startl[3];
						start*=standard.length;
						clr=(int)start;
						start=start%1;
						start*=2;
						voice=(int)start;
						start=start%1;
						start*=2;
						ssh=(int)start;
						start=start%1;
						start*=2;
						lr=(int)start;
						start=start%1;
						out=out.concat(letters[hiss[ssh][voice]]+letters[standard[clr][voice]]+letters[special[lr]]);
					}
				}
				else//two consonants
				{
					boolean done=false;
					int i=0;
					start/=startl[2];
					
					while(!done)
					{
						
					if(start<twoconsstart[i])
					{
						done=true;
						start/=twoconsstart[i];
						
						if(i!=6)
						{
							start*=standard.length;
							clr=(int)start;
							start=start%1;
						}
						start*=2;
						voice=(int)start;
						start=start%1;
						if(i==3||(i>4))
						{
							start*=2;
							ssh=(int)start;
							start=start%1;
						}
						if(i==1||i==3||i==6)
						{
							start*=2;
							lr=(int)start;
							start=start%1;
						}
						if(i<3)
							out.concat(letters[standard[clr][voice]]);
						else if(i==3)
							out=out.concat(letters[standard[clr][voice+3]]);
						else
							out=out.concat(letters[hiss[ssh][voice]]);
						//2ndletter
						if(i==4)
							out=out.concat(letters[standard[clr][voice]]);
						else if(i==5)
							out=out.concat(letters[standard[clr][2]]);
						else if(i==0)
							out=out.concat(letters[standard[clr][voice+3]]);
						else if (i==2)
							out=out.concat(letters[hiss[ssh][voice]]);
						else out=out.concat(letters[special[lr]]);
					}
					else
					{
						start-=twoconsstart[i];
						i++;
					}
					}
						
				}
			}
			else//length1 starts
			{
				start*=cons/startl[1];
				out=out.concat(letters[(int)start]);
				start=start%1;
			}
		}
		//vowel
		start*=vowel.length;
		out=out.concat(letters[vowel[(int)start]]);
		start=start%1;
		//end
		if(start>startl[0])
		{
			start-=startl[0];
			if(start>startl[1])
			{
				start-=startl[1];
				if(start>startl[2])
				{
					start-=startl[2];
					//3consonants
					{
						start/=startl[3];
						start*=standard.length;
						clr=(int)start;
						start=start%1;
						start*=2;
						voice=(int)start;
						start=start%1;
						start*=2;
						ssh=(int)start;
						start=start%1;
						start*=2;
						lr=(int)start;
						start=start%1;
						out=out.concat(letters[special[lr]]+letters[standard[clr][voice]]+letters[hiss[ssh][voice]]);
					}
				}
				else//two consonants
				{
					boolean done=false;
					int i=0;
					start/=startl[2];
					while(!done)
					{
						//System.out.println(i+","+start);
					if(start<twoconsstart[i])
					{
						done=true;
						start/=twoconsstart[i];
						
						if(i!=6)
						{
							start*=standard.length;
							clr=(int)start;
							start=start%1;
						}
						start*=2;
						voice=(int)start;
						start=start%1;
						if(i==3||(i>4))
						{
							start*=2;
							ssh=(int)start;
							start=start%1;
						}
						if(i==1||i==3||i==6)
						{
							start*=2;
							lr=(int)start;
							start=start%1;
						}
						if(i<3)
							out=out.concat(letters[standard[clr][voice]]);
						else if(i==3)
							out=out.concat(letters[standard[clr][voice+3]]);
						else
							out=out.concat(letters[hiss[ssh][voice]]);
						//2ndletter
						if(i==4)
							out=out.concat(letters[standard[clr][voice]]);
						else if(i==5)
							out=out.concat(letters[standard[clr][2]]);
						else if(i==0)
							out=out.concat(letters[standard[clr][voice+3]]);
						else if (i==2)
							out=out.concat(letters[hiss[ssh][voice]]);
						else out=out.concat(letters[special[lr]]);
					}
					else
					{
						start-=twoconsstart[i];
						i++;
					}
					}
						
				}
			}
			else//length1 starts
			{
				start*=cons/startl[1];
				out=out.concat(letters[(int)start]);
			}
		}
		return out;
	}
	void saturize() 
	{
		Color c=new Color(color);
		if (!c.equals(new Color(127,127,127)))
		{
		double red=c.getRed()-127.5,green=c.getGreen()-127.5,blue=c.getBlue()-127.5,max=Math.max(Math.abs(red), Math.max(Math.abs(green), Math.abs(blue)));
		red=red*127.5/max+128;green=green*127.5/max+128;blue=blue*127.5/max+128;
		color=new Color((int)red,(int)green,(int)blue).getRGB();
		}
	}
	@Override
	public int compareTo(java.lang.Object o) 
	{
		Particle ob=(Particle)o;
		return (int)Math.signum(mass-ob.mass);
	}

	public void print() 
	{
		System.out.print("loc{");
		for(int i=0;i<dim;i++)
			System.out.print(loc[i]+", ");
		System.out.print("}, v{");
		for(int i=0;i<dim;i++)
			System.out.print(v[i]+", ");
		System.out.println("}, m="+mass+", r="+radius);
		
	}
	//If the force is a first derivative instead of second.
	public void shove(double[] acc, double t) 
	{
		for(int i=0;i<dim;i++)
		{
			loc[i]+=t*acc[i];
			v[i]=acc[i];
		}	
		
	}
	public double torusdistance(double[]dir, double norm,double[]size)
	{
		double out=0;
		for(int i=0;i<dim;i++)
		{
			out+=Math.pow(Math.abs(dir[i]), norm);
		}
		out=Math.pow(out, 1/norm);
		
		return out;
	}
	public void setCharge(double e)
	{
		charge=e;
	}
	public double[] gravAcc(Particle obj)
	{
		double[]acc=new double[dim];
		double d=distance(obj,Cluster.metric);
		for(int k=0;k<dim;k++)
		acc[k]+=-g*obj.mass*Math.pow(d,-(rulescope+1))*(loc[k]-obj.loc[k]);
		
		return acc;
	}
	public static double[] torusgravAcc(Particle obj, double[] dir, boolean relativistic,int[] torussize, double mass2)
	{
		double[]acc=new double[dim];
		double r=(Math.pow(torussize[0], 4)*Math.pow(Math.sin(dir[0]*Math.PI/torussize[0]), 2)+
				Math.pow(torussize[1], 4)*Math.pow(Math.sin(dir[1]*Math.PI/torussize[1]), 2))*2/Math.PI;
		for(int k=0;k<dim;k++)
		{
			acc[k]=+g*obj.mass*Math.pow(torussize[k],3)*Math.sin(dir[k]/torussize[k]*Math.PI*2)/r;
			if(relativistic)
				acc[k]*=mass2;
		}
		return acc;
	}
	//For torus,klein bottle &rp2
	public static double[]coverAcc(Particle obj, double[] dir, boolean relativistic,int[] torussize, boolean[] mirror, double mass2)
	{
		double[] d=dir.clone(),
				out=new double[dim];
		int[]ts=torussize.clone();
		if(mirror[0])
		{
			ts[0]*=2;
			d[0]=d[0]-torussize[0];
			d[1]+=2*obj.loc[1];
			out=add(out,torusgravAcc(obj,d,relativistic,ts,mass2));

			if(mirror[1])
			{
				ts[1]*=2;
				d[0]+=2*obj.loc[0];
				d[1]=d[1]-torussize[1];
				out=add(out,torusgravAcc(obj,d,relativistic,ts,mass2));
			}
			d=dir.clone();
		}
		if(mirror[1])
		{
			ts[1]=torussize[1]*2;
			d[0]+=2*obj.loc[0];
			d[1]=d[1]-torussize[1];
			out=add(out,torusgravAcc(obj,d,relativistic,ts,mass2));
			d=dir.clone();
		}
		out=add(out,torusgravAcc(obj,dir,relativistic,ts,mass2));
		return out;
	}
	
	static double[] add(double[] a, double[] b) 
	{
		double[] out=new double[a.length];
		for(int i=0;i<a.length;i++)
		{
			out[i]=a[i]+b[i];
		}
		return out;
	}
	public double[] elAcc(Particle obj)
	{
		double[]acc=new double[dim];
		double d=distance(obj,Cluster.metric);
		for(int k=0;k<dim;k++)
		acc[k]+=kc*obj.charge*charge/mass*Math.pow(d,-(rulescope+1))*(loc[k]-obj.loc[k]);
		
		return acc;
	}
	
	public double[] relgravForce(Particle obj)
	{
		double[]force=new double[dim];
		double d=distance(obj,Cluster.metric);
		for(int k=0;k<dim;k++)
		force[k]+=-g*mass*obj.mass*Math.pow(d,-(rulescope+1))*(loc[k]-obj.loc[k]);
		//System.out.println("Force="+force[0]+","+force[1]);
		return force;
	}
	public void relPull(double[] force, double t) 
	{
	
		System.out.println("acc="+force[0]+","+force[1]);
		for(int i=0;i<dim;i++)
		{
			v[i]+=t*force[i];
		}
		double pnorm=vdifference(zro(),metric);
		for(int i=0;i<dim;i++)
		{
			double speed=v[i]*c/Math.pow(Math.pow(pnorm, metric)+Math.pow(mass*c, metric), 1.0/metric);
			loc[i]+=t*speed;
			//System.out.println("speed"+speed);
		}
	}
	public double[] direction(Particle particle) 
	{
		double[] out=new double[dim];
		for(int i=0;i<dim;i++)
			out[i]=particle.loc[i]-loc[i];
		return out;
	}
	public double sphericalDistance(Particle obj)
	{
		double out=Math.acos(Math.sin(loc[1])*Math.sin(obj.loc[1])*Math.cos(loc[0]-obj.loc[0])+Math.cos(loc[1])*Math.cos(obj.loc[1]));
		
		return out;
	}
	public double[] spherAcc(Particle obj, double[] direction, boolean relativistic,double r) 
	{
		double[]acc=new double[dim];
		double denom=-r*(1-Math.sin(loc[1])*Math.sin(obj.loc[1])*Math.cos(loc[0]-obj.loc[0])-Math.cos(loc[1])*Math.cos(obj.loc[1]))/4/g/obj.mass;
		acc[0]=(Math.sin(obj.loc[1])*Math.sin(loc[0]-obj.loc[0]))/denom;
		acc[1]=(Math.cos(obj.loc[1])*Math.sin(loc[1])-Math.sin(obj.loc[1])*Math.cos(loc[1])*Math.cos(loc[0]-obj.loc[0]))/denom;
		
	//	System.out.println("denom="+g);
		return acc;
	}
	
	//For e-g spherical situations
	public void adddirections(double[] acc, double[] spherAcc) 
	{
		// TODO Auto-generated method stub
		
	}
	public void spherpull(double[] acc, double t,double r) 
	{
		//System.out.print("acc=");
		//print(acc);
		double[] veucl=new double[3],leucl=new double[3];
		for(int i=0;i<2;i++)
		{
			v[i]+=t*acc[i];
		}
	
		veucl[0]=(Math.cos(loc[1])*Math.cos(loc[0])*v[1]-Math.sin(loc[0])*v[0]);
		veucl[1]=Math.sin(loc[0])*Math.cos(loc[1])*v[1]+Math.cos(loc[0])*v[0];
		veucl[2]=-Math.sin(loc[1])*v[1];
		leucl[0]=r*Math.cos(loc[0])*Math.sin(loc[1]);
		leucl[1]=r*Math.sin(loc[0])*Math.sin(loc[1]);
		leucl[2]=r*Math.cos(loc[1]);
		double vnorm=norm(veucl),gamma=Math.sqrt(vnorm*vnorm+c*c);
	
		/*System.out.print(", v_euc=");
		print(veucl);
		System.out.print(", loc_eucl=");
		print(leucl);*/
		double vangle=(vnorm*c/gamma/r-Math.PI)%(2*Math.PI)+Math.PI;
		boolean zro=(vnorm==0);

		if((!zro)&&Math.abs(vangle)!=Math.PI/2)
		{
		for(int i=0;i<3;i++)
		{
			leucl[i]+=veucl[i]*t*Math.tan(vangle)/vnorm;
		}}
		if(vangle>Math.PI/2||vangle<-Math.PI/2)
		{
			for(int i=0;i<3;i++)
			leucl[i]*=(-1);
		}
		double lnorm=norm(leucl),
				scalar=0;
		
		for(int i=0;i<3;i++)
		{
			leucl[i]*=r/lnorm;	
			scalar+=leucl[i]*veucl[i]/r;
		}
		if(vnorm!=0)
	{for(int i=0;i<3;i++)
		{
			veucl[i]-=scalar*leucl[i]/r;
		}
		double newvnorm=norm(veucl);
		for(int i=0;i<3;i++)
		{
			veucl[i]*=vnorm/newvnorm;
			if(vangle>Math.PI/2||vangle<-Math.PI/2)veucl[i]*=(-1);
		}}
		/*System.out.print(vnorm+"="+norm(veucl)+", v_euc=");
		print(veucl);
		System.out.print(", loc_eucl=");
		print(leucl);*/
		
		loc[0]=Math.atan2(leucl[1], leucl[0]);
		loc[1]=Math.acos(leucl[2]/r);
		v[0]=-Math.sin(loc[0])*veucl[0]+Math.cos(loc[0])*veucl[1];
		v[1]=Math.cos(loc[0])*Math.cos(loc[1])*veucl[0]+Math.sin(loc[0])*Math.cos(loc[1])*veucl[1]-Math.sin(loc[1])*veucl[2];
	

	}
	private static void print(double[] v) 
	{
		System.out.print("(");
		for(int i=0;i<v.length;i++)
			System.out.print(v[i]+",");
		System.out.println(")");
		
	}
	static double norm(double[] v) 
	{
		double out=0;
		for(int i=0;i<v.length;i++)
			out+=v[i]*v[i];
		out=Math.sqrt(out);
		return out;
	}
	public static Particle spherzero() 
	{
		double[] loc= {0,Math.PI/2},
				v= {0,0};
		
		return new Particle(loc,v,0,0,background);
	}

	//inner energy
	public double u()
	{
		if(spherical)
		return mass*(r*r+mass*g/2*(1+Math.cos(radius))*(Math.log(1-Math.cos(radius))-(1+Math.cos(radius))/(1-Math.cos(radius))));
		
		else return mass*(c*c);
	}
	//Just changing location by vel without affegting speed
	public void spherpush(double[] vel) 
	{
		double[]leucl=new double[3],veucl=new double[3],veucl2=new double[3];//veucl: how do I move? veucl2: what is my v-vector?
		veucl[0]=(Math.cos(loc[1])*Math.cos(loc[0])*vel[1]-Math.sin(loc[0])*vel[0]);
		veucl[1]=Math.sin(loc[0])*Math.cos(loc[1])*vel[1]+Math.cos(loc[0])*vel[0];
		veucl[2]=-Math.sin(loc[1])*vel[1];
		veucl2[0]=(Math.cos(loc[1])*Math.cos(loc[0])*v[1]-Math.sin(loc[0])*v[0]);
		veucl2[1]=Math.sin(loc[0])*Math.cos(loc[1])*v[1]+Math.cos(loc[0])*v[0];
		veucl2[2]=-Math.sin(loc[1])*v[1];
		leucl[0]=r*Math.cos(loc[0])*Math.sin(loc[1]);
		leucl[1]=r*Math.sin(loc[0])*Math.sin(loc[1]);
		leucl[2]=r*Math.cos(loc[1]);
		double vnorm=norm(veucl),
				vnorm2=norm(veucl2),
				gamma=Math.sqrt(vnorm*vnorm+c*c);
	
		/*System.out.print(", v_euc=");
		print(veucl);
		System.out.print(", loc_eucl=");
		print(leucl);*/
		double vangle=(vnorm*c/gamma/r-Math.PI)%(2*Math.PI)+Math.PI;
		boolean zro=(vnorm==0);

		if((!zro)&&Math.abs(vangle)!=Math.PI/2)
		{
		for(int i=0;i<3;i++)
		{
			leucl[i]+=veucl[i]*Math.tan(vangle)/vnorm;
		}}
		if(vangle>Math.PI/2||vangle<-Math.PI/2)
		{
			for(int i=0;i<3;i++)
			leucl[i]*=(-1);
		}
		double lnorm=norm(leucl),
				scalar=0,scalar2=0,normvecnorm;
		double[]normvec=cross(leucl,veucl),normcomp=new double[3],parcomp=new double[3];
		normvecnorm=norm(normvec);
		for(int i=0;i<3;i++)
		{
			leucl[i]/=lnorm;	
			normvec[i]/=normvecnorm;
			scalar+=normvec[i]*veucl2[i];
			scalar2+=leucl[i]*veucl2[i];
		}
		if(vnorm2!=0)
	{for(int i=0;i<3;i++)
		{
			normcomp[i]=scalar*normvec[i];
			parcomp[i]=veucl2[i]-normcomp[i];
		}
		double parnorm=norm(parcomp);
		for(int i=0;i<3;i++)
		{
			parcomp[i]-=scalar2*leucl[i];
		}
		double newparnorm=norm(parcomp);
		for(int i=0;i<3;i++)
		{
			veucl2[i]=parcomp[i]*parnorm/newparnorm;
			if(vangle>Math.PI/2||vangle<-Math.PI/2)veucl2[i]*=(-1);
			veucl2[i]+=normcomp[i];
		}}
		
		/*System.out.print(vnorm+"="+norm(veucl)+", v_euc=");
		print(veucl);
		System.out.print(", loc_eucl=");
		print(leucl);*/
		
		loc[0]=Math.atan2(leucl[1], leucl[0]);
		loc[1]=Math.acos(leucl[2]);
		v[0]=-Math.sin(loc[0])*veucl2[0]+Math.cos(loc[0])*veucl2[1];
		v[1]=Math.cos(loc[0])*Math.cos(loc[1])*veucl2[0]+Math.sin(loc[0])*Math.cos(loc[1])*veucl2[1]-Math.sin(loc[1])*veucl2[2];
	
	}
	private void normalize(double[] vec) 
	{
		double norm=norm(vec);
		for(int i=0; i<vec.length;i++)
			vec[i]/=norm;
	}
	public static double[] cross(double[] v, double[] w)
	{
		return new double[] {v[1]*w[2]-v[2]*w[1],v[2]*w[0]-v[0]*w[2],v[0]*w[1]-v[1]*w[0]};
	}
	public double[] z() 
	{
		double phi=loc[0],r=Math.tan(loc[1]/2);
		return new double[] {Math.cos(phi)*r,Math.sin(phi)*r};
	}
	
	public void addPotential(Potential p)
	{
		
		double pnorm=vdifference(zro(),metric);
		
		//	double factor=1-pnorm/Math.pow(Math.pow(pnorm, metric)+Math.pow(mass*c, metric), 1.0/metric)*0.9,
		//			cos=Math.cos(Math.PI*factor);	
		//System.out.print("factor="+factor);
		double amount=-4*g*mass/radius/radius;//(/factor);
		for(int i=(int) Math.max((loc[0]-radius)*p.res+0.99,0);i<(loc[0]+radius)*p.res;i++)
		{
			double r=Math.sqrt(radius*radius*p.res*p.res-Math.pow(loc[0]*p.res-i, 2));
			for(int j=(int)Math.max(loc[1]*p.res-r+0.99,0);j<loc[1]*p.res+r;j++)
			{
				//if(pnorm==0||(i-loc[0]*p.res)*v[0]+(j-loc[1]*p.res)*v[1]>=pnorm*Math.sqrt(Math.pow(loc[0]*p.res-i, 2)+Math.pow(loc[1]*p.res-j, 2))*cos)
				{
					p.p[i][j]+=amount;
					//System.out.print(".");
				}
			}
		}
	}
}
