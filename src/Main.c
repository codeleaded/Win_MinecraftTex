//#include "C:/Wichtig/System/Static/Library/WindowEngine1.0.h"
#include "/home/codeleaded/System/Static/Library/WindowEngine1.0.h"
#include "/home/codeleaded/System/Static/Library/Random.h"
#include "/home/codeleaded/System/Static/Library/PerlinNoise.h"
#include "/home/codeleaded/System/Static/Library/RayCast.h"

#include "./Math3D.h"
#include "./Hitbox3D.h"


typedef struct Camera {
	mat4x4 matProj;
	vec3d vCamera;
	vec3d vVelocity;
	vec3d vLookDir;
	float fYaw;	
	float fPitch;
	vec3d vLength;
} Camera;

Camera Camera_New(vec3d p){
	Camera c;
	c.matProj = Matrix_MakeProjection(90.0f,(float)GetHeight() / (float)GetWidth(),0.1f,1000.0f);
	c.vCamera = p;
	c.vVelocity = (vec3d){ 0.0f,0.0f,0.0f,1.0f };
	c.vLookDir = (vec3d){ 0.0f,0.0f,0.0f,1.0f };
	c.fYaw = 0.0f;
	c.fPitch = 0.0f;
	c.vLength = (vec3d){ 0.0f,0.0f,0.0f,1.0f };
	return c;
}


mesh meshCube;
mat4x4 matProj;
vec3d vCamera = { 10.0f,50.0f,10.0f,1.0f };
vec3d vVelocity = { 0.0f,0.0f,0.0f,1.0f };
vec3d vLookDir;
float fYaw;	
float fPitch;	
float fTheta;

vec3d light_direction;

Sprite sprTex1;
float *pDepthBuffer;

vec3d vLength = { 0.5f,1.8f,0.5f,1.0f };
Vector Cubes;
Vec2 MouseBefore = { 0.0f,0.0f };

mat4x4 matView;

char OnGround = 0;
int Mode = 0;
int Menu = 0;

typedef unsigned char Block;

#define BLOCK_ERROR	255

#define BLOCK_VOID	0
#define BLOCK_DIRT	1
#define BLOCK_GRAS	2
#define BLOCK_LEAF	3
#define BLOCK_LOG	4
#define BLOCK_SAND	5
#define BLOCK_WATER	6

#define ATLAS_VOID		255
#define ATLAS_DIRT		0
#define ATLAS_GRAS_S	1
#define ATLAS_GRAS_T	2
#define ATLAS_SAND		3
#define ATLAS_WATER		4
#define ATLAS_LOG_S		5
#define ATLAS_LOG_T		6
#define ATLAS_LEAF		7

#define WORLD_DX	50
#define WORLD_DY	40
#define WORLD_DZ	50

Block World[WORLD_DX * WORLD_DY * WORLD_DZ];

Block World_GetX(Block* World,float x,float y,float z){
	if(x<0.0f || x>=WORLD_DX) return BLOCK_ERROR;
	if(y<0.0f || y>=WORLD_DY) return BLOCK_ERROR;
	if(z<0.0f || z>=WORLD_DZ) return BLOCK_ERROR;
	return World[(int)x + (int)y * WORLD_DX + (int)z * WORLD_DX * WORLD_DY];
}
void World_SetX(Block* World,float x,float y,float z,Block b){
	if(x<0.0f || x>=WORLD_DX) return;
	if(y<0.0f || y>=WORLD_DY) return;
	if(z<0.0f || z>=WORLD_DZ) return;
	World[(int)x + (int)y * WORLD_DX + (int)z * WORLD_DX * WORLD_DY] = b;
}
Block World_Get(Block* World,vec3d p){
	return World_GetX(World,p.x,p.y,p.z);
}
void World_Set(Block* World,vec3d p,Block b){
	World_SetX(World,p.x,p.y,p.z,b);
}

int World_Height(Block* World,float x,float z){
	int h = WORLD_DY - 1;
	for(;h>=0;h--){
		if(World_GetX(World,x,h,z)!=BLOCK_VOID){
			return h;
		}
	}
	return h;
}

#define CUBE_SIDE_SOUTH		0
#define CUBE_SIDE_EAST		1
#define CUBE_SIDE_NORTH		2
#define CUBE_SIDE_WEST		3
#define CUBE_SIDE_TOP		4
#define CUBE_SIDE_BOTTOM	5


vec3d Neighbour_Side(int s){
	switch (s){
	case CUBE_SIDE_SOUTH: 	return (vec3d){ 0.0f, 0.0f,-1.0f,1.0f };
	case CUBE_SIDE_EAST: 	return (vec3d){ 1.0f, 0.0f, 0.0f,1.0f };
	case CUBE_SIDE_NORTH: 	return (vec3d){ 0.0f, 0.0f, 1.0f,1.0f };
	case CUBE_SIDE_WEST: 	return (vec3d){-1.0f, 0.0f, 0.0f,1.0f };
	case CUBE_SIDE_TOP: 		return (vec3d){ 0.0f, 1.0f, 0.0f,1.0f };
	case CUBE_SIDE_BOTTOM: 	return (vec3d){ 0.0f,-1.0f, 0.0f,1.0f };
	}
	return (vec3d){ 0.0f,0.0f,0.0f,1.0f };
}
unsigned char Block_Id(Block b,int s){
	switch (b){
	case BLOCK_VOID: 	return ATLAS_VOID;
	case BLOCK_DIRT: 	return ATLAS_DIRT;
	case BLOCK_GRAS:	return s==CUBE_SIDE_TOP ? ATLAS_GRAS_T : (s==CUBE_SIDE_BOTTOM ? ATLAS_DIRT : ATLAS_GRAS_S);
	
	case BLOCK_SAND:	return ATLAS_SAND;
	case BLOCK_WATER:	return ATLAS_WATER;
	case BLOCK_LOG: 	return s==CUBE_SIDE_TOP || s==CUBE_SIDE_BOTTOM ? ATLAS_LOG_T : ATLAS_LOG_S;
	case BLOCK_LEAF: 	return ATLAS_LEAF;
	}
	return ATLAS_VOID;
}

void World_Tree(Block* World,int x,int y,int z){
	World_SetX(World,x,y,z,BLOCK_LOG);
	World_SetX(World,x,y+1,z,BLOCK_LOG);
	World_SetX(World,x,y+2,z,BLOCK_LOG);

	World_SetX(World,x,y+3,z,BLOCK_LEAF);

	World_SetX(World,x+1,y+3,z,BLOCK_LEAF);
	World_SetX(World,x+1,y+3,z+1,BLOCK_LEAF);
	World_SetX(World,x+1,y+3,z-1,BLOCK_LEAF);
	World_SetX(World,x-1,y+3,z,BLOCK_LEAF);
	World_SetX(World,x-1,y+3,z+1,BLOCK_LEAF);
	World_SetX(World,x-1,y+3,z-1,BLOCK_LEAF);
	World_SetX(World,x,y+3,z+1,BLOCK_LEAF);
	World_SetX(World,x,y+3,z-1,BLOCK_LEAF);

	World_SetX(World,x,y+4,z,BLOCK_LEAF);
	World_SetX(World,x+1,y+4,z,BLOCK_LEAF);
	World_SetX(World,x-1,y+4,z,BLOCK_LEAF);
	World_SetX(World,x,y+4,z+1,BLOCK_LEAF);
	World_SetX(World,x,y+4,z-1,BLOCK_LEAF);
}
void World_Generate(Block* World,int dx,int dy,int dz){
	memset(World,0,dx * dy * dz);

	float Seed[dx * dz];
	for(int i = 0;i<dx * dz;i++){
		Seed[i] = (float)Random_f64_New();
	}

	float Out[dx * dz];
	PerlinNoise_2D_Buffer(dx,dz,Seed,15,0.9f,Out);

	for(int i = 0;i<dz;i++){
		for(int j = 0;j<dx;j++){
			int h = 0;
			for(;h<Out[j + i * dx] * ((float)dy * 0.8f);h++){
				World_SetX(World,j,h,i,BLOCK_DIRT);
			}
			World_SetX(World,j,h,i,BLOCK_GRAS);
		}
	}


	for(int i = 2;i<dz-2;i++){
		for(int j = 2;j<dx-2;j++){
			int h = World_Height(World,j,i) + 1;
			
			if(Random_i32_MinMax(0,40)==0){
				World_Tree(World,j,h,i);
			}
		}
	}
}

void Triangle_CalcNorm(triangle* t){
	vec3d normal, line1, line2;

	line1 = vec3d_Sub(t->p[1],t->p[0]);
	line2 = vec3d_Sub(t->p[2],t->p[0]);

	normal = vec3d_CrossProduct(line1,line2);

	t->n = vec3d_Normalise(normal);
}

void Cube_Set(triangle* trisout,vec3d p,vec3d d,Block id){
	triangle tris[12] = {
	// SOUTH
	{ 0.0f,0.0f,0.0f,1.0f,  0.0f,1.0f,0.0f,1.0f,  1.0f,1.0f,0.0f,1.0f,  0.0f,1.0f,1.0f,  0.0f,0.0f,1.0f,  1.0f,0.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_SOUTH)},
	{ 0.0f,0.0f,0.0f,1.0f,  1.0f,1.0f,0.0f,1.0f,  1.0f,0.0f,0.0f,1.0f,  0.0f,1.0f,1.0f,  1.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_SOUTH)},
	// EAST                                                      
	{ 1.0f,0.0f,0.0f,1.0f,  1.0f,1.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,1.0f,  0.0f,1.0f,1.0f,  0.0f,0.0f,1.0f,  1.0f,0.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_EAST)},
	{ 1.0f,0.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,1.0f,  1.0f,0.0f,1.0f,1.0f,  0.0f,1.0f,1.0f,  1.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_EAST)},
	// NORTH                                                     
	{ 1.0f,0.0f,1.0f,1.0f,  1.0f,1.0f,1.0f,1.0f,  0.0f,1.0f,1.0f,1.0f,  0.0f,1.0f,1.0f,  0.0f,0.0f,1.0f,  1.0f,0.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_NORTH)},
	{ 1.0f,0.0f,1.0f,1.0f,  0.0f,1.0f,1.0f,1.0f,  0.0f,0.0f,1.0f,1.0f,  0.0f,1.0f,1.0f,  1.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_NORTH)},
	// WEST                                                      
	{ 0.0f,0.0f,1.0f,1.0f,  0.0f,1.0f,1.0f,1.0f,  0.0f,1.0f,0.0f,1.0f,  0.0f,1.0f,1.0f,  0.0f,0.0f,1.0f,  1.0f,0.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_WEST)},
	{ 0.0f,0.0f,1.0f,1.0f,  0.0f,1.0f,0.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  0.0f,1.0f,1.0f,  1.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_WEST)},
	// TOP                                                       
	{ 0.0f,1.0f,0.0f,1.0f,  0.0f,1.0f,1.0f,1.0f,  1.0f,1.0f,1.0f,1.0f,  0.0f,1.0f,1.0f,  0.0f,0.0f,1.0f,  1.0f,0.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_TOP)},
	{ 0.0f,1.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,1.0f,  1.0f,1.0f,0.0f,1.0f,  0.0f,1.0f,1.0f,  1.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_TOP)},
	// BOTTOM                                                    
	{ 1.0f,0.0f,1.0f,1.0f,  0.0f,0.0f,1.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  0.0f,1.0f,1.0f,  0.0f,0.0f,1.0f,  1.0f,0.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_BOTTOM)},
	{ 1.0f,0.0f,1.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  1.0f,0.0f,0.0f,1.0f,  0.0f,1.0f,1.0f,  1.0f,0.0f,1.0f,  1.0f,1.0f,1.0f,  0.0f,0.0f,0.0f,1.0f,  WHITE,Block_Id(id,CUBE_SIDE_BOTTOM)},
	};

	for(int i = 0;i<12;i++){
		for(int j = 0;j<3;j++){
			tris[i].p[j] = vec3d_Add(p,vec3d_Make(tris[i].p[j].x * d.x,tris[i].p[j].y * d.y,tris[i].p[j].z * d.z));
		}
		Triangle_CalcNorm(&tris[i]);

		float dp = fmaxf(0.1f, vec3d_DotProduct(light_direction,tris[i].n));
		tris[i].c = Pixel_Mul(tris[i].c,Pixel_toRGBA(dp,dp,dp,1.0f));
		
		trisout[i] = tris[i];
	}
}
void MakeCube(vec3d p,vec3d d,Block id){
	triangle tris[12];
	Cube_Set(tris,p,d,id);

	for(int i = 0;i<12;i++){
		Vector_Push(&meshCube.tris,&tris[i]);
	}
}
void MakePlane(vec3d p,vec3d d,int Plane,Block id){
	triangle tris[12];
	Cube_Set(tris,p,d,id);

	for(int i = Plane*2;i<(Plane+1)*2;i++){
		Vector_Push(&meshCube.tris,&tris[i]);
	}
}
void BuildCube(vec3d p,vec3d d,Block id){
	Vector_Push(&Cubes,(Rect3[]){ { p,d } });
	MakeCube(p,d,id);
}

void Mesh_Reload(){
	Vector_Clear(&meshCube.tris);

	for(int i = 0;i<WORLD_DZ;i++){
		for(int j = 0;j<WORLD_DY;j++){
			for(int k = 0;k<WORLD_DX;k++){
				Block b = World_GetX(World,k,j,i);

				if(b!=BLOCK_VOID){
					for(int s = 0;s<6;s++){
						vec3d p = { k,j,i,1.0f };
						vec3d n = vec3d_Add(p,Neighbour_Side(s));
						
						if(World_Get(World,n)==BLOCK_VOID){
							MakePlane(p,(vec3d){ 1.0f,1.0f,1.0f },s,b);
						}
					}
				}
			}
		}
	}
}

int Cubes_Compare(const void* e1,const void* e2) {
	Rect3 r1 = *(Rect3*)e1;
	Rect3 r2 = *(Rect3*)e2;
	
	vec3d pos = vec3d_Add(vCamera,(vec3d){ vLength.x * 0.5f,vLength.y * 0.9f,vLength.z * 0.5f });
	vec3d d1 = vec3d_Sub(r1.p,pos);
    vec3d d2 = vec3d_Sub(r2.p,pos);
	return vec3d_Length(d1) == vec3d_Length(d2) ? 0 : (vec3d_Length(d1) < vec3d_Length(d2) ? 1 : -1);
}
void Cubes_Reload(Block* World){
	Vector_Clear(&Cubes);

	vec3d f = { (int)vCamera.x,(int)vCamera.y,(int)vCamera.z };
	for(int i = -2;i<2;i++){
		for(int j = -2;j<2;j++){
			for(int k = -2;k<2;k++){
				vec3d n = { k,j,i };
				vec3d r = vec3d_Add(f,n);

				Block b = World_Get(World,r);

				if(b!=BLOCK_VOID && b!=BLOCK_ERROR){
					Vector_Push(&Cubes,(Rect3[]){ { r,(vec3d){ 1.0f,1.0f,1.0f } } });
				}
			}
		}
	}

	qsort(Cubes.Memory,Cubes.size,Cubes.ELEMENT_SIZE,Cubes_Compare);
}

void World_Edit(Block* World,vec3d p,Block b){
	World_SetX(World,p.x,p.y,p.z,b);
	Mesh_Reload();
}

void int_swap(int* a,int* b){
	int c = *a;
	*a = *b;
	*b = c;
}
void float_swap(float* a,float* b){
	float c = *a;
	*a = *b;
	*b = c;
}
float float_max(float a,float b){
	return a>b ? a : b;
}

void TexturedTriangle(	int x1, int y1, float u1, float v1, float w1,
						int x2, int y2, float u2, float v2, float w2,
						int x3, int y3, float u3, float v3, float w3,
                    	Pixel col,unsigned char id,Sprite *tex
){
	if (y2 < y1){
		int_swap(&y1,&y2);
		int_swap(&x1,&x2);
		float_swap(&u1,&u2);
		float_swap(&v1,&v2);
		float_swap(&w1,&w2);
	}
	if (y3 < y1){
		int_swap(&y1,&y3);
		int_swap(&x1,&x3);
		float_swap(&u1,&u3);
		float_swap(&v1,&v3);
		float_swap(&w1,&w3);
	}
	if (y3 < y2){
		int_swap(&y2,&y3);
		int_swap(&x2,&x3);
		float_swap(&u2,&u3);
		float_swap(&v2,&v3);
		float_swap(&w2,&w3);
	}

	int dy1 = y2 - y1;
	int dx1 = x2 - x1;
	float dv1 = v2 - v1;
	float du1 = u2 - u1;
	float dw1 = w2 - w1;
	int dy2 = y3 - y1;
	int dx2 = x3 - x1;
	float dv2 = v3 - v1;
	float du2 = u3 - u1;
	float dw2 = w3 - w1;
	float tex_u, tex_v, tex_w;
	float dax_step = 0, dbx_step = 0,
		du1_step = 0, dv1_step = 0,
		du2_step = 0, dv2_step = 0,
		dw1_step=0, dw2_step=0;
	if (dy1) dax_step = dx1 / (float)fabsf(dy1);
	if (dy2) dbx_step = dx2 / (float)fabsf(dy2);
	if (dy1) du1_step = du1 / (float)fabsf(dy1);
	if (dy1) dv1_step = dv1 / (float)fabsf(dy1);
	if (dy1) dw1_step = dw1 / (float)fabsf(dy1);
	if (dy2) du2_step = du2 / (float)fabsf(dy2);
	if (dy2) dv2_step = dv2 / (float)fabsf(dy2);
	if (dy2) dw2_step = dw2 / (float)fabsf(dy2);
	
	if (dy1){
		for (int i = y1; i <= y2; i++){
			int ax = x1 + (float)(i - y1) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;
			float tex_su = u1 + (float)(i - y1) * du1_step;
			float tex_sv = v1 + (float)(i - y1) * dv1_step;
			float tex_sw = w1 + (float)(i - y1) * dw1_step;
			float tex_eu = u1 + (float)(i - y1) * du2_step;
			float tex_ev = v1 + (float)(i - y1) * dv2_step;
			float tex_ew = w1 + (float)(i - y1) * dw2_step;
			
			if (ax > bx){
				int_swap(&ax,&bx);
				float_swap(&tex_su,&tex_eu);
				float_swap(&tex_sv,&tex_ev);
				float_swap(&tex_sw,&tex_ew);
			}
			
			tex_u = tex_su;
			tex_v = tex_sv;
			tex_w = tex_sw;
			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;
			for (int j = ax; j < bx; j++){
				tex_u = (1.0f - t) * tex_su + t * tex_eu;
				tex_v = (1.0f - t) * tex_sv + t * tex_ev;
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
				if (tex_w > pDepthBuffer[i*GetWidth() + j]){
					int subx = (int)id % 16;
					int suby = (int)id / 16;
					Pixel sample = Sprite_SampleSub(tex,tex_u / tex_w, tex_v / tex_w,subx,suby,16,16);
					//sample = Pixel_Mul(col,sample);
					Draw(j,i,sample);
					pDepthBuffer[i*GetWidth() + j] = tex_w;
				}
				t += tstep;
			}
		}
	}

	dy1 = y3 - y2;
	dx1 = x3 - x2;
	dv1 = v3 - v2;
	du1 = u3 - u2;
	dw1 = w3 - w2;
	if (dy1) dax_step = dx1 / (float)fabsf(dy1);
	if (dy2) dbx_step = dx2 / (float)fabsf(dy2);
	du1_step = 0, dv1_step = 0;
	if (dy1) du1_step = du1 / (float)fabsf(dy1);
	if (dy1) dv1_step = dv1 / (float)fabsf(dy1);
	if (dy1) dw1_step = dw1 / (float)fabsf(dy1);
	
	if(dy1){
		for(int i = y2; i <= y3; i++){
			int ax = x2 + (float)(i - y2) * dax_step;
			int bx = x1 + (float)(i - y1) * dbx_step;
			float tex_su = u2 + (float)(i - y2) * du1_step;
			float tex_sv = v2 + (float)(i - y2) * dv1_step;
			float tex_sw = w2 + (float)(i - y2) * dw1_step;
			float tex_eu = u1 + (float)(i - y1) * du2_step;
			float tex_ev = v1 + (float)(i - y1) * dv2_step;
			float tex_ew = w1 + (float)(i - y1) * dw2_step;
			
			if(ax > bx){
				int_swap(&ax,&bx);
				float_swap(&tex_su,&tex_eu);
				float_swap(&tex_sv,&tex_ev);
				float_swap(&tex_sw,&tex_ew);
			}

			tex_u = tex_su;
			tex_v = tex_sv;
			tex_w = tex_sw;
			float tstep = 1.0f / ((float)(bx - ax));
			float t = 0.0f;
			for(int j = ax; j < bx; j++){
				tex_u = (1.0f - t) * tex_su + t * tex_eu;
				tex_v = (1.0f - t) * tex_sv + t * tex_ev;
				tex_w = (1.0f - t) * tex_sw + t * tex_ew;
				if (tex_w > pDepthBuffer[i*GetWidth() + j])
				{
					int subx = (int)id % 16;
					int suby = (int)id / 16;
					Pixel sample = Sprite_SampleSub(tex,tex_u / tex_w, tex_v / tex_w,subx,suby,16,16);
					//sample = Pixel_Mul(col,sample);
					Draw(j,i,sample);
					pDepthBuffer[i * GetWidth() + j] = tex_w;
				}
				t += tstep;
			}
		}	
	}		
}

void Triangles_Project(){
	Vector vecTrianglesToRaster = Vector_New(sizeof(triangle));
	for(int i = 0;i<meshCube.tris.size;i++){
		triangle tri = *(triangle*)Vector_Get(&meshCube.tris,i);
		
		vec3d vCameraRay = vec3d_Sub(tri.p[0],vCamera);
		if (vec3d_DotProduct(tri.n,vCameraRay) < 0.0f){
			tri.p[0] = Matrix_MultiplyVector(matView,tri.p[0]);
			tri.p[1] = Matrix_MultiplyVector(matView,tri.p[1]);
			tri.p[2] = Matrix_MultiplyVector(matView,tri.p[2]);
			
			int nClippedTriangles = 0;
			triangle clipped[2];
			nClippedTriangles = Triangle_ClipAgainstPlane(vec3d_Make(0.0f,0.0f,0.1f),vec3d_Make(0.0f,0.0f,1.0f),&tri,&clipped[0],&clipped[1]);
			
			for (int n = 0; n < nClippedTriangles; n++){
				triangle triProjected;
				triProjected.p[0] = Matrix_MultiplyVector(matProj,clipped[n].p[0]);
				triProjected.p[1] = Matrix_MultiplyVector(matProj,clipped[n].p[1]);
				triProjected.p[2] = Matrix_MultiplyVector(matProj,clipped[n].p[2]);
				triProjected.c = clipped[n].c;
				triProjected.id = clipped[n].id;
				triProjected.t[0] = clipped[n].t[0];
				triProjected.t[1] = clipped[n].t[1];
				triProjected.t[2] = clipped[n].t[2];

				triProjected.t[0].u = triProjected.t[0].u / triProjected.p[0].w;
				triProjected.t[1].u = triProjected.t[1].u / triProjected.p[1].w;
				triProjected.t[2].u = triProjected.t[2].u / triProjected.p[2].w;
				triProjected.t[0].v = triProjected.t[0].v / triProjected.p[0].w;
				triProjected.t[1].v = triProjected.t[1].v / triProjected.p[1].w;
				triProjected.t[2].v = triProjected.t[2].v / triProjected.p[2].w;
				triProjected.t[0].w = 1.0f / triProjected.p[0].w;
				triProjected.t[1].w = 1.0f / triProjected.p[1].w;
				triProjected.t[2].w = 1.0f / triProjected.p[2].w;
				
				triProjected.p[0] = vec3d_Div(triProjected.p[0],triProjected.p[0].w);
				triProjected.p[1] = vec3d_Div(triProjected.p[1],triProjected.p[1].w);
				triProjected.p[2] = vec3d_Div(triProjected.p[2],triProjected.p[2].w);
				
				triProjected.p[0].x *= -1.0f;
				triProjected.p[1].x *= -1.0f;
				triProjected.p[2].x *= -1.0f;
				triProjected.p[0].y *= -1.0f;
				triProjected.p[1].y *= -1.0f;
				triProjected.p[2].y *= -1.0f;
				
				vec3d vOffsetView = vec3d_Make(1.0f,1.0f,0.0f);
				triProjected.p[0] = vec3d_Add(triProjected.p[0],vOffsetView);
				triProjected.p[1] = vec3d_Add(triProjected.p[1],vOffsetView);
				triProjected.p[2] = vec3d_Add(triProjected.p[2],vOffsetView);
				triProjected.p[0].x *= 0.5f * (float)GetWidth();
				triProjected.p[0].y *= 0.5f * (float)GetHeight();
				triProjected.p[1].x *= 0.5f * (float)GetWidth();
				triProjected.p[1].y *= 0.5f * (float)GetHeight();
				triProjected.p[2].x *= 0.5f * (float)GetWidth();
				triProjected.p[2].y *= 0.5f * (float)GetHeight();

				Vector_Push(&vecTrianglesToRaster,&triProjected);
			}			
		}
	}
	
	Clear(BLACK);
	free(pDepthBuffer);
	pDepthBuffer = malloc(sizeof(float) * GetWidth() * GetHeight());
	memset(pDepthBuffer,0,sizeof(float) * GetWidth() * GetHeight());
		
	for(int i = 0;i<vecTrianglesToRaster.size;i++){
		triangle triToRaster = *(triangle*)Vector_Get(&vecTrianglesToRaster,i);
		
		triangle clipped[2];
		Vector listTriangles = Vector_New(sizeof(triangle));
		
		Vector_Push(&listTriangles,&triToRaster);
		int nNewTriangles = 1;
		for (int p = 0; p < 4; p++){
			int nTrisToAdd = 0;
			while (nNewTriangles > 0){
				triangle test = *(triangle*)Vector_Get(&listTriangles,0);
				Vector_Remove(&listTriangles,0);
				nNewTriangles--;
				
				switch (p)
				{
				case 0:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3d){ 0.0f, 0.0f, 0.0f }, 					(vec3d){ 0.0f, 1.0f, 0.0f },&test,&clipped[0],&clipped[1]); break;
				case 1:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3d){ 0.0f, (float)GetHeight() - 1, 0.0f }, 	(vec3d){ 0.0f,-1.0f, 0.0f },&test,&clipped[0],&clipped[1]); break;
				case 2:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3d){ 0.0f, 0.0f, 0.0f }, 					(vec3d){ 1.0f, 0.0f, 0.0f },&test,&clipped[0],&clipped[1]); break;
				case 3:	nTrisToAdd = Triangle_ClipAgainstPlane((vec3d){ (float)GetWidth() - 1, 0.0f, 0.0f }, 	(vec3d){-1.0f, 0.0f, 0.0f },&test,&clipped[0],&clipped[1]); break;
				}
				
				for (int w = 0; w < nTrisToAdd; w++)
					Vector_Push(&listTriangles,&clipped[w]);
			}
			nNewTriangles = listTriangles.size;
		}
		
		for (int i = 0;i<listTriangles.size;i++){
			triangle t = *(triangle*)Vector_Get(&listTriangles,i);

			if(Mode==0){
				TexturedTriangle(
					t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
					t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
					t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w, 
					t.c,t.id,&sprTex1
				);
			}
			if(Mode==1) RenderTriangleWire(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),t.c,1.0f);
			if(Mode==2){
				RenderTriangle(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),t.c);
			}
			if(Mode==3){
				TexturedTriangle(
					t.p[0].x, t.p[0].y, t.t[0].u, t.t[0].v, t.t[0].w,
					t.p[1].x, t.p[1].y, t.t[1].u, t.t[1].v, t.t[1].w,
					t.p[2].x, t.p[2].y, t.t[2].u, t.t[2].v, t.t[2].w, 
					t.c,t.id,&sprTex1
				);
				RenderTriangleWire(((Vec2){ t.p[0].x, t.p[0].y }),((Vec2){ t.p[1].x, t.p[1].y }),((Vec2){ t.p[2].x, t.p[2].y }),WHITE,1.0f);
			}
			
			//RenderTriangle(((Vec2){t.p[0].x,t.p[0].y}),((Vec2){t.p[1].x,t.p[1].y}),((Vec2){t.p[2].x,t.p[2].y}),t.c);
			//RenderTriangleWire(((Vec2){t.p[0].x,t.p[0].y}),((Vec2){t.p[1].x,t.p[1].y}),((Vec2){t.p[2].x,t.p[2].y}),t.c,1.0f);
			//RenderTriangleWire(((Vec2){t.p[0].x,t.p[0].y}),((Vec2){t.p[1].x,t.p[1].y}),((Vec2){t.p[2].x,t.p[2].y}),WHITE,1.0f);
		}
		Vector_Free(&listTriangles);
	}
	Vector_Free(&vecTrianglesToRaster);
}

void Stand(vec3d* Data){
	Data->y = 0.0f;
	OnGround = 1;
}

char World_Void(Block* World,Vec3 p){
	Block b = World_GetX(World,p.x,p.y,p.z);
	return b==BLOCK_VOID || b==BLOCK_ERROR;
}

void Menu_Set(int m){
	if(Menu==0 && m==1){
		AlxWindow_Mouse_SetInvisible(&window);
		SetMouse((Vec2){ GetWidth() / 2,GetHeight() / 2 });
	}
	if(Menu==1 && m==0){
		AlxWindow_Mouse_SetVisible(&window);
	}
	
	MouseBefore = GetMouse();
	Menu = m;
}

void Setup(AlxWindow* w){
	Menu_Set(1);
	RGA_Get(6969);

	meshCube = (mesh){ Vector_New(sizeof(triangle)) };
	Cubes = Vector_New(sizeof(Rect3));

	light_direction = vec3d_Normalise(vec3d_Make(0.4f,0.5f,-0.6f));

	pDepthBuffer = malloc(sizeof(float) * GetWidth() * GetHeight());
	
	sprTex1 = Sprite_Load("./data/Atlas.png");

	matProj = Matrix_MakeProjection(90.0f, (float)GetHeight() / (float)GetWidth(), 0.1f, 1000.0f);

	World_Generate(World,WORLD_DX,WORLD_DY,WORLD_DZ);

	Mesh_Reload();
}

void Update(AlxWindow* w){
	//w->ElapsedTime = 0.05;
	
	/*if(Stroke(ALX_MOUSE_L).PRESSED){
		MouseBefore = GetMouse();
	}
	if(Stroke(ALX_MOUSE_L).DOWN){
		if(GetMouse().x!=MouseBefore.x || GetMouse().y!=MouseBefore.y){
			Vec2 d = Vec2_Sub(GetMouse(),MouseBefore);
			Vec2 a = Vec2_Mulf(Vec2_Div(d,(Vec2){ window.Width,window.Height }),2 * PI);
	
			fYaw += a.x;
			fPitch += a.y;
	
			MouseBefore = GetMouse();
		}
	}*/
	if(Menu==1){
		if(GetMouse().x!=MouseBefore.x || GetMouse().y!=MouseBefore.y){
			Vec2 d = Vec2_Sub(GetMouse(),MouseBefore);
			Vec2 a = Vec2_Mulf(Vec2_Div(d,(Vec2){ window.Width,window.Height }),2 * F32_PI);
	
			fYaw += a.x;
			fPitch += a.y;
	
			SetMouse((Vec2){ GetWidth() / 2,GetHeight() / 2 });
			MouseBefore = GetMouse();
		}
	}
	
	if(Stroke(ALX_KEY_ESC).PRESSED)
		Menu_Set(!Menu);

	if(Stroke(ALX_KEY_UP).DOWN)
		fPitch -= 2.0f * w->ElapsedTime;

	if(Stroke(ALX_KEY_DOWN).DOWN)
		fPitch += 2.0f * w->ElapsedTime;

	if(Stroke(ALX_KEY_LEFT).DOWN)
		fYaw -= 2.0f * w->ElapsedTime;

	if(Stroke(ALX_KEY_RIGHT).DOWN)
		fYaw += 2.0f * w->ElapsedTime;

	if(Stroke(ALX_KEY_Z).PRESSED)
		Mode = Mode < 3 ? Mode+1 : 0;

	if(Stroke(ALX_KEY_R).DOWN)
		//if(OnGround) 
			vVelocity.y = 5.0f;
	
	//if(Stroke(ALX_KEY_F).DOWN)
	//	vVelocity.y = -5.0f;

	//if(Stroke(ALX_KEY_R).RELEASED || Stroke(ALX_KEY_F).RELEASED)
	//	vVelocity.y = 0.0f;

	mat4x4 matCameraRot = Matrix_MakeRotationY(fYaw);
	vec3d vForward = Matrix_MultiplyVector(matCameraRot,vec3d_Make(0.0f,0.0f,1.0f));
	vec3d vLeft = vec3d_PerpY(vForward);
	
	if(Stroke(ALX_KEY_W).DOWN){
		vVelocity.x += vForward.x * 20.0f * w->ElapsedTime;
		vVelocity.z += vForward.z * 20.0f * w->ElapsedTime;
	}
	if(Stroke(ALX_KEY_S).DOWN){
		vVelocity.x -= vForward.x * 20.0f * w->ElapsedTime;
		vVelocity.z -= vForward.z * 20.0f * w->ElapsedTime;
	}
	if(Stroke(ALX_KEY_A).DOWN){
		vVelocity.x -= vLeft.x * 20.0f * w->ElapsedTime;
		vVelocity.z -= vLeft.z * 20.0f * w->ElapsedTime;
	}
	if (Stroke(ALX_KEY_D).DOWN){
		vVelocity.x += vLeft.x * 20.0f * w->ElapsedTime;
		vVelocity.z += vLeft.z * 20.0f * w->ElapsedTime;
	}

	Vec2 v = { vVelocity.x,vVelocity.z };
	Vec2 d = Vec2_Norm(v);

	float drag = OnGround ? 12.0f : 10.0f;
	Vec2 da = Vec2_Norm(v);
	v = Vec2_Sub(v,Vec2_Mulf(d,drag * w->ElapsedTime));

	if(F32_Sign(v.x)!=F32_Sign(da.x) || F32_Sign(v.y)!=F32_Sign(da.y)){
		v.x = 0.0f;
		v.y = 0.0f;
	}

	float maxVelo = OnGround ? 4.0f : 6.0f;
	if(Vec2_Mag(v)>maxVelo){
		v = Vec2_Mulf(d,maxVelo);
	}

	vVelocity.x = v.x;
	vVelocity.z = v.y;

	vVelocity = vec3d_Add(vVelocity,vec3d_Mul((vec3d){ 0.0f,-10.0f,0.0f,1.0f },w->ElapsedTime));
	vCamera = vec3d_Add(vCamera,vec3d_Mul(vVelocity,w->ElapsedTime));

	Cubes_Reload(World);
	OnGround = 0;
	for(int i = 0;i<Cubes.size;i++){
		vec3d pos = { vLength.x * 0.5f,vLength.y * 0.9f,vLength.z * 0.5f };

		Rect3 r1 = *(Rect3*)Vector_Get(&Cubes,i);
		Rect3 r2 = (Rect3){ vec3d_Sub(vCamera,pos),vLength };
		Rect3_Static(&r2,r1,&vVelocity,(void (*[])(void*)){ NULL,NULL,NULL,NULL,NULL,(void*)Stand });
		vCamera = vec3d_Add(r2.p,pos);
	}

	float Border = F32_PI * 0.5f - 0.00001;
	if(fPitch<-Border) fPitch = -Border;
	if(fPitch>Border) fPitch = Border;

	vec3d vUp = vec3d_Make(0.0f,1.0f,0.0f);
	vec3d vTarget = vec3d_Make(0.0f,0.0f,1.0f);
	mat4x4 matCameraRotX = Matrix_MakeRotationX(fPitch);
	vLookDir = Matrix_MultiplyVector(matCameraRotX,vTarget);
	vLookDir = Matrix_MultiplyVector(matCameraRot,vLookDir);
	
	vTarget = vec3d_Add(vCamera, vLookDir);
	mat4x4 matCamera = Matrix_PointAt(vCamera, vTarget, vUp);
	matView = Matrix_QuickInverse(matCamera);

	if(Stroke(ALX_MOUSE_L).PRESSED){
		Vec3 c = (Vec3){ vCamera.x,vCamera.y,vCamera.z };
		RayCast_TileMap(World,(void*)World_Void,c,(Vec3){ vLookDir.x,vLookDir.y,vLookDir.z },0.01f,4.0f,&c);
		if(c.x!=vCamera.x || c.y!=vCamera.y || c.z!=vCamera.z)
			World_Edit(World,(vec3d){ c.x,c.y,c.z },BLOCK_VOID);
	}
	if(Stroke(ALX_MOUSE_R).PRESSED){
		Vec3 c = (Vec3){ vCamera.x,vCamera.y,vCamera.z };
		RayCast_TileMap_N(World,(void*)World_Void,c,(Vec3){ vLookDir.x,vLookDir.y,vLookDir.z },0.01f,4.0f,&c);
		if(c.x!=vCamera.x || c.y!=vCamera.y || c.z!=vCamera.z)
			World_Edit(World,(vec3d){ c.x,c.y,c.z },BLOCK_DIRT);
	}

	Clear(LIGHT_BLUE);
	Triangles_Project();

	String str = String_Format("X: %f, Y: %f, Z: %f, Size: %d",vCamera.x,vCamera.y,vCamera.z,meshCube.tris.size);
	RenderCStrSize(str.Memory,str.size,0,0,RED);
	String_Free(&str);
}

void Delete(AlxWindow* w){
	free(pDepthBuffer);
	Sprite_Free(&sprTex1);

	Vector_Free(&meshCube.tris);
	Vector_Free(&Cubes);

	AlxWindow_Mouse_SetVisible(&window);
}

int main(){
    if(Create("3D Test no Tex",2500,1200,1,1,Setup,Update,Delete))
        Start();
    return 0;
}