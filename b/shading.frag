#version 450

layout(location = 0) in vec2 inTexcoord; 
layout(location =0) out vec4 fragColor;


layout(std140, binding = 0) uniform UBO
{
    vec2 viewportSize;
    vec2 mousePos;
    int mouseButtonAction[3];
    float cameraDist;
};

uniform float focusStrength = 0.1; // 控制散焦强度，可以根据需求调整
vec4 color = vec4(1.0);
vec4 outColor = vec4(1.0);

// 计算散焦效果
vec4 applyDof(vec4 color, float depth)
{
    float focusDist = cameraDist;
    float focusRange = 2.0; // 控制焦距范围，可以根据需求调整

    float dofAmount = clamp(abs(depth - focusDist) / focusRange, 0.0, 1.0);
    dofAmount = pow(dofAmount, focusStrength);

    return mix(color, vec4(0.0, 0.0, 0.0, 1.0), dofAmount);
}

#define PI 3.14159265358979323846
#define CAMERA_DIST 3

float maxcomp(vec2 v){
	return max(v.x,v.y);
}

float maxcomp(vec3 v){
	return max(max(v.x,v.y),v.z);
}
float maxcomp(vec4 v){
	return max(max(max(v.x,v.y),v.z),v.w);
}
float mincomp(vec2 v){
	return min(v.x,v.y);
}
float mincomp(vec3 v){
	return min(min(v.x,v.y),v.z);
}
float mincomp(vec4 v){
	return min(min(min(v.x,v.y),v.z),v.w);
}

float hash11( float n )    // in [0,1]
{
    return fract(sin(n*13.)*43758.5453);
}

vec2 hash22( vec2 p ) {
    return fract(sin(vec2(dot(p,vec2(127.1,311.7)),dot(p,vec2(269.5,183.3))))*43758.5453);
}


#define DEG2RAD(a) a*PI/180.

mat3 identity(inout mat3 m){
    m[0]=vec3(0.);
    m[1]=vec3(0.);
    m[2]=vec3(0.);
    m[0][0]=1.;
    m[1][1]=1.;
    m[2][2]=1.;
    return m;
}
//rotate with x y z
mat3 rotate(float angle, vec3 v){
    float a = DEG2RAD(angle);
	float c = cos(a);
	float s = sin(a);

    vec3 axis=normalize(v);
    vec3 temp=(1. - c) * axis;

    mat3 Rotate;
    Rotate[0][0] = c + temp[0] * axis[0];
	Rotate[0][1] = temp[0] * axis[1] + s * axis[2];
	Rotate[0][2] = temp[0] * axis[2] - s * axis[1];

	Rotate[1][0] = temp[1] * axis[0] - s * axis[2];
	Rotate[1][1] = c + temp[1] * axis[1];
	Rotate[1][2] = temp[1] * axis[2] + s * axis[0];

	Rotate[2][0] = temp[2] * axis[0] + s * axis[1];
	Rotate[2][1] = temp[2] * axis[1] - s * axis[0];
	Rotate[2][2] = c + temp[2] * axis[2];
	return Rotate;
}


mat3 scale(vec3 factor){
    mat3 m;
    identity(m);
    m[0][0]=factor[0];
    m[1][1]=factor[1];
    m[2][2]=factor[2];
    return m;
}


//eye matrix
mat3 lookAt(vec3 pos,vec3 lookat,vec3 up){
    mat3 m;
    identity(m);
    
    //eye matrix
    vec3 ww=-normalize(lookat-pos);
    vec3 uu=normalize(cross(up,ww));
    vec3 vv=normalize(cross(ww,uu));
    m[0]=uu;
    m[1]=vv;
    m[2]=ww;
    return m;    
}

mat3 surfaceTBN(vec3 normal){
    mat3 m;
    identity(m);
    vec3 randomVec = normal==vec3(0,1,0)?normalize(vec3(.234234,1.13565,5.42342)):vec3(0,1,0);
    // create TBN change-of-basis matrix: from tangent-space to view-space
    vec3 tangent = normalize(randomVec - normal * dot(randomVec, normal));
    vec3 bitangent = cross(normal, tangent);
    mat3 TBN = mat3(tangent, bitangent, normal);
    return TBN;
}

/////////////////////////////////////////////////////////////
#define RAY_MAXDEPTH 100.
#define RAY_INTERVAL 0.001
#define ZERO (min(iFrame,0))


#define ROOM 0
#define SQUARE 1
#define SPHERE1 2
#define CUBE1 3
#define LIGHT 4

//refractive index
#define RF_AIR 1.0 
#define RF_WATER 1.33 
#define RF_GLASS 1.5

const float sphereR=3.;
const int LIGHTSNUM=1;

const mat4 identityM=mat4(vec4(1.,0.,0.,0.),vec4(0.,1.,0.,0.),vec4(0.,0.,1.,0.),vec4(0.,0.,0.,1.));
const vec3 triangluarLightVertices[]=vec3[3](vec3(-3.,19.9,0.),vec3(3.,19.9,0.),vec3(0.,19.9,3.));

struct Material{
    vec3 albedo;
	vec3 kd;	//diffuse
	vec3 ks;	//specular
    float shinness;	//shinees
};


struct Light{
	vec3 pos;
    vec3 intensity;
};

struct Obj{
	float depth;
    int id;
    vec3 pos;
    vec3 N;
};

struct Para{
    Obj obj;
    vec3 V;
    vec3 k; //light intensity coefficient
};
const vec3 a_glass=vec3(.3,0.2,0.3); //attenuation constant

    
Obj intersectDepth(in Obj A, in Obj B) {
    if(A.depth>B.depth)return A;
    return B;
}

Obj unionDepth(Obj A, Obj B) {
    if(A.depth<B.depth)return A;
    return B;
}

Material getMaterialById(int id){
    switch(id){
        case ROOM: //house
        return Material(vec3(1.),vec3(0.6),vec3(0.4),1000.);
        case SPHERE1: //sphere
        return Material(vec3(1.,0.,0.),vec3(0.6),vec3(0.4),100.);
        case CUBE1: 
        return Material(vec3(1.,1.,0.),vec3(0.6),vec3(0.9),100.);
        case SQUARE:
        return Material(vec3(1.,0.,0.),vec3(0.6),vec3(0.4),100.);   
    }
    //no intersection
    return Material(vec3(0.),vec3(0.),vec3(0.),0.);
}
float checkersMod( in vec3 p )
{
    vec3 q = floor(p);
    return clamp(mod(q.x+q.y+q.z,2.),0.,1.);
}

float intersectCube(vec3 ro,vec3 rd){
    //cube size s
    vec3 temp[2];
    temp[0]=(-1.-ro)/rd;
    temp[1]=(1.-ro)/rd;
    vec3 tmin= min(temp[0],temp[1]);
    vec3 tmax= max(temp[0],temp[1]);

    float tMaxMin=maxcomp(tmin);
    float tMinMax=mincomp(tmax);
    
    if(tMinMax>=tMaxMin&&tMinMax>RAY_INTERVAL)return (tMaxMin>RAY_INTERVAL?tMaxMin:tMinMax);
    return 10000.;
}
Obj intersectCube1(vec3 ro,vec3 rd){
    Obj res;
    vec3 offset=vec3(0,1,3);
    vec3 scaleFactor=vec3(1);
    ro-=offset;

    mat3 M=rotate(30.,vec3(0.,1.,0.))*scale(scaleFactor);
    mat3 Mi=inverse(M);
    ro=Mi*ro;
    rd=Mi*rd;

    float t=intersectCube(ro,rd);
    vec3 pos=ro+t*rd;
    res.N=normalize(M*round(pos*.50001));
    res.pos=M*pos+offset;
    res.id=CUBE1;
    res.depth=t;
    return res;
}

float intersectSquare(vec3 ro,vec3 rd){
    float t=-ro.y/rd.y;
    vec3 pos=ro+t*rd;
    if(maxcomp(abs(pos.xz)-vec2(1))<0.&&t>RAY_INTERVAL)return t;
    return 10000.;
}
float intersectSphere(in vec3 ro,in vec3 rd){
    float h;
    //ray intersection
    vec3 ce=ro;
    float a=dot(rd,rd);
    float b=dot(ce,rd);
    float c=dot(ce,ce)-1.;
    h=b*b-a*c;
    if(h>0.){
        float sqrth=sqrt(h);
    	float h1=(-b-sqrth)/a;
        float h2=(-b+sqrth)/a;
        if(h2>=h1&&h2>RAY_INTERVAL)return (h1>RAY_INTERVAL?h1:h2);
    }
    return 10000.;
}
Obj intersectSphere1(vec3 ro,vec3 rd){
    Obj res;
    vec3 offset=vec3(0.,1.1,0.);
    ro-=offset;
    float t=intersectSphere(ro,rd);
    vec3 pos=ro+t*rd;
    res.N=pos;
    res.pos=pos+offset;
    res.id=SPHERE1;
    res.depth=t;
    return res;
}
Obj intersectSquare1(vec3 ro,vec3 rd){
    Obj res;
    vec3 offset=vec3(0.,0.1,0.);
    vec3 scaleFactor=vec3(7.);
    ro-=offset;

    mat3 M=scale(scaleFactor);
    mat3 Mi=inverse(M);
    ro=Mi*ro;
    rd=Mi*rd;

    float t=intersectSquare(ro,rd);
    vec3 pos=ro+t*rd;
    res.N=vec3(0.,1.,0.);
    res.pos=M*pos+offset;
    res.id=SQUARE;
    res.depth=t;
    return res;

}
//setup scene
Obj intersectScene(vec3 ro,vec3 rd){
    Obj res=Obj(10000.,-1,vec3(0.),vec3(0.));
    res=unionDepth(res,intersectSquare1(ro,rd));
    res=unionDepth(res,intersectSphere1(ro,rd));
    res=unionDepth(res,intersectCube1(ro,rd));
    return res;
}
float calcShadow(vec3 ro, vec3 rd,float tmin,float tmax,float k){
	float res=1.;
    Obj temp;
    temp=intersectScene(ro,rd);
    if(temp.depth>tmin&&temp.depth<tmax)return 0.;
    
    return 1.;
}

float SchlickApproximation(float cosTheta){
    const float R0=0.0016;// glass
    // return R0+(1.-R0)*pow(1.-cosTheta,5.);
    return R0 + (1.0 - R0) * exp2((-5.55473 * cosTheta - 6.98316) * cosTheta);
}

vec3 calcShadow(vec3 ro,vec3 N,vec3 L,float tmin,float tmax){
    const int sampleNum=16;

	vec3 res=vec3(0.);
    vec3 rd[sampleNum];
    mat3 TBN=surfaceTBN(L);
    
    for(int i=0;i<sampleNum;++i){
        vec2 randomVec=hash22(ro.xy+vec2(hash11(float(i)),hash11(13.+float(i)*17.)));

		vec3 tempVec=vec3(0.,0.,1.);
        tempVec.xy=(randomVec-0.5)*0.3;
        rd[i]=normalize(L+TBN*tempVec);
        //rd[i]=L;
        float curRF=RF_AIR,nextRF=RF_AIR;

        vec3 k=vec3(1.);
        Para curPara;
        curPara.obj=intersectScene(ro,rd[i]);
        curPara.V=rd[i];
        curPara.k=k;

        for(int i=0;i<3;++i){
            if(curRF==RF_GLASS)curPara.obj.N=-curPara.obj.N;
            if(curPara.obj.depth>tmin&&curPara.obj.depth<tmax){
                //hit object
                if(curPara.obj.id==SPHERE1||curPara.obj.id==CUBE1){
                    nextRF= curRF==RF_AIR?RF_GLASS:RF_AIR;    //hit glass
                    vec3 t=refract(curPara.V,curPara.obj.N,curRF/nextRF);
                    Obj nextObj=intersectScene(curPara.obj.pos,t); //refraction
                    float cosTheta=-dot(curPara.V,curPara.obj.N);
                    float R=SchlickApproximation(cosTheta);
                    vec3 kR;
                    if(any(notEqual(t,vec3(0.)))){
                        curPara.V=t;
                        //into glass
                        if(nextRF==RF_GLASS)k=exp(-a_glass*nextObj.depth);
                        curRF=nextRF;
                        kR=curPara.k*k*(1.-R);
                        curPara.k=kR;
                        if(nextObj.depth>tmin&&nextObj.depth<tmax){
                            curPara.obj=nextObj;
                        }else{break;} //attenuation to 0 or hit nothing
                    }else{
                        
                        k=vec3(1.0);
                        //total reflection or too small
                        vec3 r=reflect(curPara.V,curPara.obj.N);
                        curPara.V=r;
                        nextObj=intersectScene(curPara.obj.pos,r); //reflection
                        if(nextRF==RF_AIR)k=exp(-a_glass*nextObj.depth);
                        kR=curPara.k*k*R;
                        if(any(greaterThan(kR,vec3(0.01)))&&nextObj.depth>tmin&&nextObj.depth<tmax){
                            curPara.obj=nextObj;
                        }else {break;}
                    }; //total reflection or too small
                }else {curPara.k=vec3(0.);break;};    //not glass
            }
        }
		float cTheta=max(dot(curPara.V,L),0.);
        cTheta=max((cTheta-0.9)*10.,0.);
        res+=cTheta*curPara.k;
    }
    res*=2./float(sampleNum);

    return res;
}

//with shadow
vec3 BlinnPhongShading(in Light lights[LIGHTSNUM],in Material material,vec3 pos,vec3 V,vec3 N){
    //shadow
	vec3 col=vec3(0.);
    for(int i=0;i<LIGHTSNUM;++i){
        #if 0
        vec3 tempL=lights[i].pos-pos;
        vec3 L=normalize(tempL);
        vec3 intensity=lights[i].intensity/pow(length(tempL),2.);

        #else
        vec3 L=normalize(lights[i].pos);
        vec3 intensity=lights[i].intensity; 	//no attenuation
        #endif
        vec3 H=normalize(L-V);
        //float k_shadow=calcShadow(pos,L,RAY_INTERVAL+0.001,RAY_MAXDEPTH,8.);
        vec3 k_shadow=calcShadow(pos,N,L,RAY_INTERVAL,RAY_MAXDEPTH);
        col+=k_shadow*intensity*(material.albedo*material.kd*max(0.,dot(N,L))+material.ks*pow(max(0.,dot(N,H)),material.shinness));    
    }
    return (col);
}

void setCamera(vec2 uv,out vec3 ro,out vec3 rd){
    uv-=0.5;
    if(viewportSize.x>viewportSize.y) 
        uv.x*=float(viewportSize.x)/viewportSize.y;
    else
        uv.y*=viewportSize.y/viewportSize.x;

    vec2 mo=mousePos.xy+0.5;
	float theta=PI*(mo.y-0.5);
	float phi=2.*PI*mo.x-PI/2.;
	//eye parameters
	float ans=0;//.*iTime;

	ro=cameraDist*vec3(cos(phi)*cos(theta),0.2,sin(phi)*cos(theta));
	vec3 up=vec3(0,1,0);
	vec3 center=vec3(0,0,0);	//look at point
	
	//eye matrix
	mat3 V=lookAt(ro,center,up);
	//cast a view ray
	rd=normalize(V*vec3(uv,-1.));	//world space
}

vec3 render(vec3 ro,vec3 rd){

 	vec3 bgColor=vec3(0.5);
    
    Light light1[LIGHTSNUM];
    light1[0]=Light(vec3(1.,0.5,1.),vec3(1.));
    vec3 ambient=vec3(0.1);

    vec3 col=vec3(0.);
    float curRF=RF_AIR,nextRF=RF_AIR;
    
    Para interObjs[10];
    vec3 k=vec3(1.);
    interObjs[0].k=vec3(1.);

    
    interObjs[0].obj=intersectScene(ro,rd); 
    if(interObjs[0].obj.depth>0.0001&&interObjs[0].obj.depth<RAY_MAXDEPTH){
        interObjs[0].V=rd;
        for(int i=0;i<10&&i>-1;){
            Para curPara=interObjs[i];
            if(curRF==RF_GLASS)curPara.obj.N=-curPara.obj.N;
            Material curM=getMaterialById(curPara.obj.id);
            vec3 k=vec3(1.);
            //col=m.albedo;
            if(curPara.obj.id==LIGHT){--i;}
            else if(curPara.obj.id==SPHERE1||curPara.obj.id==CUBE1){
                //col=curPara.obj.N;
                //break;
                nextRF= curRF==RF_AIR?RF_GLASS:RF_AIR;    //hit glass
                //glass
                vec3 t=refract(curPara.V,curPara.obj.N,curRF/nextRF);
                Obj nextObj=intersectScene(curPara.obj.pos,t); //refraction
                float cosTheta=-dot(curPara.V,curPara.obj.N);
                float R=SchlickApproximation(cosTheta);
                vec3 kR;
                if(any(notEqual(t,vec3(0.)))&&cosTheta>.005){
                    //into
                    if(nextRF==RF_GLASS)k=exp(-a_glass*nextObj.depth);
                     curRF=nextRF;
                    //push refraction
                    kR=curPara.k*k*(1.-R);
                    if(nextObj.depth>0.001&&nextObj.depth<RAY_MAXDEPTH){
                        if(any(greaterThan(kR,vec3(0.005)))){
                            interObjs[i].k=kR;
                            interObjs[i].obj=nextObj;
                            interObjs[i].V=t;
                        }else --i;
                    }else{ 
						 col+=kR*bgColor;--i;}
                }else {--i;R=1.;}
                //reflection
                k=vec3(1.);
                kR=curPara.k*R;
                // 
                if(any(greaterThan(kR,vec3(.001)))&&cosTheta>.005){
                    //col=vec3(1.0);
                    vec3 r=reflect(curPara.V,curPara.obj.N);
                    nextObj=intersectScene(curPara.obj.pos,r); //reflection
                    if(nextRF==RF_AIR)k=exp(-a_glass*nextObj.depth);
                    if(nextObj.depth>0.0001&&nextObj.depth<RAY_MAXDEPTH){
                        ++i;
                        interObjs[i].k=kR*k;
                        interObjs[i].obj=nextObj;
                        interObjs[i].V=r;
                    }else if(nextObj.depth>=RAY_MAXDEPTH){
						//bgColor=pow(bgColor,vec3(2.2));
                        col+=kR*k*bgColor;
                    }
                };
            }
            else{
                if(curPara.obj.id==SQUARE){
                    curM.albedo=vec3(checkersMod(curPara.obj.pos));
                }
                vec3 color=BlinnPhongShading(light1,curM,curPara.obj.pos,curPara.V,curPara.obj.N);
                color+=curM.albedo*ambient;
                col+=interObjs[i].k*color;
                --i;
            }
        }
    }else{ 
		//bgColor=pow(bgColor,vec3(2.2));    
        col=bgColor;
    }
    return col;
}

void main() {
    vec3 ro = vec3(0.0, 0.0, -CAMERA_DIST);
    vec2 uv = (inTexcoord * 2.0 - 1.0) * vec2(viewportSize.x / viewportSize.y, 1.0);
    vec3 rd = normalize(vec3(uv, 1.0));
    Obj hitObj = intersectScene(ro, rd);

    if (hitObj.id != -1) {
        Material mat = getMaterialById(hitObj.id);
        vec3 color = mat.albedo * hitObj.N;
        fragColor = applyDof(vec4(color, 1.0), hitObj.depth);
    } else {
        fragColor = vec4(0.0, 0.0, 0.0, 1.0);
    }
}