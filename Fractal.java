import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import javax.imageio.ImageIO;
//**********************************************************************************************
//Author: Non-Euclidean Dreamer
//3-Body-Problem Fractal
//**********************************************************************************************



public class Fractal 
{
	static String name="Fractal";
	static int orange=new Color(255,128,0).getRGB(),
			green=Color.green.getRGB(),
			black=Color.black.getRGB();
	public static int c=black,a=0, Tmax=8000;
	public static double[][]loc= {{-10,0},{10,0}},v= {{0,0.2},{0,-0.2}};
	public static double[]mass= {1,1},radius= {1,1},
			l= {0,a/100.0},
			focus= {0,0};
	public static int[] color= {orange,green};
	public static double m=0.1,r=1,
			scale=0.25, T=0.01;
	
	
	public static int width=540, height=540;
	
	public static void main(String[] args) 
	{	if(args.length>0)
		name=args[0];
		for(int i=1;i<args.length;i++)
		{
			name.concat(args[i]);
			String string=args[i];
			if(string.substring(0,4).equals("loc0"))loc[0][Integer.parseInt(string.substring(4,5))]=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,4).equals("loc1"))loc[1][Integer.parseInt(string.substring(4,5))]=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,4).equals("loc2"))l[Integer.parseInt(string.substring(4,5))]=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,4).equals("vel0"))v[0][Integer.parseInt(string.substring(4,5))]=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,4).equals("vel1"))v[1][Integer.parseInt(string.substring(4,5))]=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,5).equals("mass0"))mass[0]=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,5).equals("mass1"))mass[1]=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,5).equals("mass2"))m=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,5).equals("radi0"))radius[0]=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,5).equals("radi1"))radius[1]=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,5).equals("radi2")) { r=Double.parseDouble(string.substring(5,string.length()));System.out.println("r"+r);}
			else if(string.substring(0,5).equals("focus"))focus[Integer.parseInt(string.substring(5,6))]=Double.parseDouble(string.substring(6,string.length()));
			else if(string.substring(0,5).equals("scale"))scale=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,5).equals("tstep"))T=Double.parseDouble(string.substring(5,string.length()));
			else if(string.substring(0,4).equals("tmax")) {Tmax=Integer.parseInt(string.substring(4,string.length()));System.out.println("tmax="+Tmax);}
			else{System.out.println("Failed parsing "+string);}
		}
	//	Scanner in=new Scanner(System.in);
	//	System.out.println("Choose parameter value:");
		//int a=in.nextInt();System.out.println("Value="+a);
		
		Particle.g=1;
		BufferedImage canvas=new BufferedImage(width,height,BufferedImage.TYPE_4BYTE_ABGR);;
		/*try {
			canvas = ImageIO.read(new File("Fractal0.png"));System.out.println("success!");
		} catch (IOException e1) {
			System.out.println("Error");
			e1.printStackTrace();
		}*/
		for(int i=0;i<width;i++)
		{System.out.println(i+"-th row");
			for(int j=0;j<height;j++)
			{
				
				ArrayList<Particle>objects=new ArrayList<Particle>();
								objects.add(new Particle(l,new double[] {(2*i-width)*scale/width+focus[0],(2*j-height)*scale/width+focus[1]},m,r,c));

				for(int k=0;k<loc.length;k++) 
					objects.add(new Particle(loc[k],v[k],mass[k],radius[k],color[k]));
				Cluster cluster=new Cluster(objects);
				//cluster.print(); 
				boolean notdone=true;int col,out=black;int t=0;
				for(t=0; t<Tmax;t++)
				{	notdone=false;
					col=cluster.update(T);
				
					//break-off condition for escaping objects(Think through...)

				
					if(col!=0) {out=Cluster.blur(black,col,(t/1.01)/Tmax);/*if(t!=0)System.out.println(t);*/t=Tmax;}
					//else if(Math.abs(objects.get(0).loc[0])+Math.abs(objects.get(0).loc[1])>500)t=Tmax;
				}//cluster.print();
				canvas.setRGB(i,j,out);
			}
		File outputfile = new File(name+a+".png");
		try 
		{
			ImageIO.write(canvas, "png", outputfile);
		} catch (IOException e) {
			System.out.println("IOException");System.out.println("couldn't print");
			e.printStackTrace(); 
		}}
	}

}
