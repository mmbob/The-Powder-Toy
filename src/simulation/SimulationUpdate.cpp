#include <cmath>
#if !defined(_MSC_VER)
#include <strings.h>
#else
#include <windows.h>
#endif
#include "Config.h"
#include "Simulation.h"
#include "Elements.h"
#include "Air.h"
#include "Gravity.h"
#include "elements/Element.h"

#include "cat/LuaScriptInterface.h"
#include "cat/LuaScriptHelper.h"

#include <SDL_thread.h>

#include <ppl.h>

float heat_map[YRES][XRES];

SDL_mutex* updateMutex = nullptr;

int part_temp_transition(Simulation* sim, int i, int x, int y);
int part_pressure_transition(Simulation* sim, int i, int x, int y);
void part_movement(Simulation* sim, int i, int x, int y, int nt, int surround_space);
float part_conduct_heat(Simulation* sim, int i, int x, int y, int surround[8]);
bool get_gravity_at(Simulation* sim, int x, int y, float elementGravity, float* gravX, float* gravY);

void Simulation::update_particles_i(int start, int inc)
{
	if (updateMutex == nullptr)
		updateMutex = SDL_CreateMutex();

	auto safe_kill_part = [=](int i)
	{
		SDL_LockMutex(updateMutex);
		kill_part(i);
		SDL_UnlockMutex(updateMutex);
	};

	currentTick++;

	bool lighting_ok = true;

	if (lighting_recreate > 0)
	{
		for (int i = 0; i<=parts_lastActiveIndex; i++)
		{
			if (parts[i].type==PT_LIGH && parts[i].tmp2>0)
			{
				lighting_ok = false;
				break;
			}
		}
	}

	if (lighting_ok)
		lighting_recreate--;

	if (lighting_recreate < 0)
		lighting_recreate=1;

	if (lighting_recreate > 21)
		lighting_recreate=21;

	//if (sys_pause&&!framerender)//do nothing if paused
	//	return;

	if (force_stacking_check || (rand() % 10)==0)
	{
		force_stacking_check = 0;
		bool excessive_stacking_found = false;
		for (int y = 0; y < YRES; y++)
		{
			for (int x = 0; x < XRES; x++)
			{
				// Use a threshold, since some particle stacking can be normal (e.g. BIZR + FILT)
				// Setting pmap_count[y][x] > NPART means BHOL will form in that spot
				if (pmap_count[y][x]>5)
				{
					if (bmap[y/CELL][x/CELL]==WL_EHOLE)
					{
						// Allow more stacking in E-hole
						if (pmap_count[y][x]>1500)
						{
							pmap_count[y][x] = pmap_count[y][x] + NPART;
							excessive_stacking_found = true;
						}
					}
					else if (pmap_count[y][x]>1500 || (rand()%1600)<=(pmap_count[y][x]+100))
					{
						pmap_count[y][x] = pmap_count[y][x] + NPART;
						excessive_stacking_found = true;
					}
				}
			}
		}
		if (excessive_stacking_found)
		{
			for (int i = 0; i <= parts_lastActiveIndex; i++)
			{
				if (parts[i].type)
				{
					int t = parts[i].type;
					int x = (int)(parts[i].x+0.5f);
					int y = (int)(parts[i].y+0.5f);
					if (x>=0 && y>=0 && x<XRES && y<YRES && !(elements[t].Properties&TYPE_ENERGY))
					{
						if (pmap_count[y][x]>=NPART)
						{
							if (pmap_count[y][x]>NPART)
							{
								create_part(i, x, y, PT_NBHL);
								parts[i].temp = MAX_TEMP;
								parts[i].tmp = pmap_count[y][x]-NPART;//strength of grav field
								if (parts[i].tmp>51200) parts[i].tmp = 51200;
								pmap_count[y][x] = NPART;
							}
							else
							{
								SDL_LockMutex(updateMutex);
								kill_part(i);
								SDL_UnlockMutex(updateMutex);
							}
						}
					}
				}
			}
		}
	}

	if (elementCount[PT_LOVE] > 0 || elementCount[PT_LOLZ] > 0) //LOVE and LOLZ element handling
	{
		for (int ny = 0; ny<YRES-4; ny++)
		{
			for (int nx = 0; nx<XRES-4; nx++)
			{
				int r = pmap[ny][nx];
				if (!r)
				{
					continue;
				}
				else if ((ny<9||nx<9||ny>YRES-7||nx>XRES-10)&&(parts[r>>8].type==PT_LOVE||parts[r>>8].type==PT_LOLZ))
					kill_part(r>>8);
				else if (parts[r>>8].type==PT_LOVE)
				{
					Element_LOVE::love[nx/9][ny/9] = 1;
				}
				else if (parts[r>>8].type==PT_LOLZ)
				{
					Element_LOLZ::lolz[nx/9][ny/9] = 1;
				}
			}
		}
		for (int nx=9; nx<=XRES-18; nx++)
		{
			for (int ny=9; ny<=YRES-7; ny++)
			{
				if (Element_LOVE::love[nx/9][ny/9]==1)
				{
					for ( int nnx=0; nnx<9; nnx++)
						for ( int nny=0; nny<9; nny++)
						{
							if (ny+nny>0&&ny+nny<YRES&&nx+nnx>=0&&nx+nnx<XRES)
							{
								int rt=pmap[ny+nny][nx+nnx];
								if (!rt&&Element_LOVE::RuleTable[nnx][nny]==1)
									create_part(-1,nx+nnx,ny+nny,PT_LOVE);
								else if (!rt)
									continue;
								else if (parts[rt>>8].type==PT_LOVE&&Element_LOVE::RuleTable[nnx][nny]==0)
									kill_part(rt>>8);
							}
						}
				}
				Element_LOVE::love[nx/9][ny/9]=0;
				if (Element_LOLZ::lolz[nx/9][ny/9]==1)
				{
					for ( int nnx=0; nnx<9; nnx++)
						for ( int nny=0; nny<9; nny++)
						{
							if (ny+nny>0&&ny+nny<YRES&&nx+nnx>=0&&nx+nnx<XRES)
							{
								int rt=pmap[ny+nny][nx+nnx];
								if (!rt&&Element_LOLZ::RuleTable[nny][nnx]==1)
									create_part(-1,nx+nnx,ny+nny,PT_LOLZ);
								else if (!rt)
									continue;
								else if (parts[rt>>8].type==PT_LOLZ&&Element_LOLZ::RuleTable[nny][nnx]==0)
									kill_part(rt>>8);

							}
						}
				}
				Element_LOLZ::lolz[nx/9][ny/9]=0;
			}
		}
	}

	//wire!
	if(elementCount[PT_WIRE] > 0)
	{
		for (int nx=0; nx<XRES; nx++)
		{
			for (int ny=0; ny<YRES; ny++)
			{
				int r = pmap[ny][nx];
				if (!r)
					continue;
				if(parts[r>>8].type==PT_WIRE)
					parts[r>>8].tmp=parts[r>>8].ctype;
			}
		}
	}

	if (Element_PPIP::ppip_changed)
	{
		for (int i=0; i<=parts_lastActiveIndex; i++)
		{
			if (parts[i].type==PT_PPIP)
			{
				parts[i].tmp |= (parts[i].tmp&0xE0000000)>>3;
				parts[i].tmp &= ~0xE0000000;
			}
		}
		Element_PPIP::ppip_changed = 0;
	}

	//game of life!
	if (elementCount[PT_LIFE]>0&&++CGOL>=GSPEED)//GSPEED is frames per generation
	{
		CGOL=0;
		//TODO: maybe this should only loop through active particles
		for (int ny=CELL; ny<YRES-CELL; ny++)
		{//go through every particle and set neighbor map
			for (int nx=CELL; nx<XRES-CELL; nx++)
			{
				int r = pmap[ny][nx];
				if (!r)
				{
					gol[ny][nx] = 0;
					continue;
				}
				if ((r&0xFF)==PT_LIFE)
				{
					int golnum = parts[r>>8].ctype+1;
					if (golnum<=0 || golnum>NGOL) {
						kill_part(r>>8);
						continue;
					}
					gol[ny][nx] = golnum;
					if (parts[r>>8].tmp == grule[golnum][9]-1) {
						for ( int nnx=-1; nnx<2; nnx++)
						{
							for ( int nny=-1; nny<2; nny++)//it will count itself as its own neighbor, which is needed, but will have 1 extra for delete check
							{
								int adx = ((nx+nnx+XRES-3*CELL)%(XRES-2*CELL))+CELL;
								int ady = ((ny+nny+YRES-3*CELL)%(YRES-2*CELL))+CELL;
								int rt = pmap[ady][adx];
								if (!rt || (rt&0xFF)==PT_LIFE)
								{
									//the total neighbor count is in 0
									gol2[ady][adx][0] ++;
									//insert golnum into neighbor table
									for ( int i=1; i<9; i++)
									{
										if (!gol2[ady][adx][i])
										{
											gol2[ady][adx][i] = (golnum<<4)+1;
											break;
										}
										else if((gol2[ady][adx][i]>>4)==golnum)
										{
											gol2[ady][adx][i]++;
											break;
										}
									}
								}
							}
						}
					} else {
						parts[r>>8].tmp --;
					}
				}
			}
		}
		for (int ny=CELL; ny<YRES-CELL; ny++)
		{ //go through every particle again, but check neighbor map, then update particles
			for (int nx=CELL; nx<XRES-CELL; nx++)
			{
				int r = pmap[ny][nx];
				if (r && (r&0xFF)!=PT_LIFE)
					continue;
				int neighbors = gol2[ny][nx][0];
				if (neighbors)
				{
					int golnum = gol[ny][nx];
					if (!r)
					{
						//Find which type we can try and create
						int creategol = 0xFF;
						for ( int i=1; i<9; i++)
						{
							if (!gol2[ny][nx][i]) break;
							golnum = (gol2[ny][nx][i]>>4);
							if (grule[golnum][neighbors]>=2 && (gol2[ny][nx][i]&0xF)>=(neighbors%2)+neighbors/2)
							{
								if (golnum<creategol) creategol=golnum;
							}
						}
						if (creategol<0xFF)
							create_part(-1, nx, ny, PT_LIFE|((creategol-1)<<8));
					}
					else if (grule[golnum][neighbors-1]==0 || grule[golnum][neighbors-1]==2)//subtract 1 because it counted itself
					{
						if (parts[r>>8].tmp==grule[golnum][9]-1)
							parts[r>>8].tmp --;
					}
					for ( int z = 0; z<9; z++)
						gol2[ny][nx][z] = 0;//this improves performance A LOT compared to the memset, i was getting ~23 more fps with this.
				}
				//we still need to kill things with 0 neighbors (higher state life)
				if (r && parts[r>>8].tmp<=0)
						kill_part(r>>8);
			}
		}
		//memset(gol2, 0, sizeof(gol2));
	}
	if (ISWIRE>0)//wifi channel reseting
	{
		for ( int q = 0; q<(int)(MAX_TEMP-73.15f)/100+2; q++)
		{
			wireless[q][0] = wireless[q][1];
			wireless[q][1] = 0;
		}
		ISWIRE--;
	}

	elementRecount |= !(currentTick%180);
	if(elementRecount)
	{
		std::fill(elementCount, elementCount+PT_NUM, 0);
	}

	concurrency::parallel_for(0, parts_lastActiveIndex + 1, [&](int i)
	{
		Particle* part = &parts[i];
		int t = part->type;
		if (t == 0)
			return;
		Element* element = &elements[t];
		if (t < 0 || t >= PT_NUM || !element->Enabled)
		{
			safe_kill_part(i);
			return;
		}

		int x = (int) (part->x + 0.5f);
		int y = (int) (part->y + 0.5f);

		int elem_properties = element->Properties;
		if (part->life > 0 && (elem_properties & PROP_LIFE_DEC))
		{
			// automatically decrease life
			part->life--;
			if (part->life <= 0 && (elem_properties & (PROP_LIFE_KILL_DEC | PROP_LIFE_KILL)))
			{
				// kill on change to no life
				safe_kill_part(i);
				return;
			}
		}
		else if (part->life <= 0 && (elem_properties & PROP_LIFE_KILL))
		{
			// kill if no life
			safe_kill_part(i);
			return;
		}
			
		int wallType = bmap[y / CELL][y / CELL];
		if (x < CELL || y < CELL || x >= XRES-CELL || y >= YRES-CELL ||
			(wallType &&
			(wallType==WL_WALL ||
			wallType==WL_WALLELEC ||
			wallType==WL_ALLOWAIR ||
			(wallType==WL_DESTROYALL) ||
			(wallType==WL_ALLOWLIQUID && element->Falldown!=2) ||
			(wallType==WL_ALLOWSOLID && element->Falldown!=1) ||
			(wallType==WL_ALLOWGAS && !(element->Properties&TYPE_GAS)) || //&& elements[t].Falldown!=0 && parts[i].type!=PT_FIRE && parts[i].type!=PT_SMKE && parts[i].type!=PT_CFLM) ||
			(wallType==WL_ALLOWENERGY && !(element->Properties&TYPE_ENERGY)) ||
			(wallType==WL_DETECT && (t==PT_METL || t==PT_SPRK)) ||
			(wallType==WL_EWALL && !emap[y/CELL][x/CELL])) && (t!=PT_STKM) && (t!=PT_STKM2) && (t!=PT_FIGH)))
		{
			safe_kill_part(i);
			return;
		}

		int surround[8];
		int surround_index = 0;
		for (int nx=-1; nx<2; nx++)
		{
			for (int ny=-1; ny<2; ny++)
			{
				if (nx||ny)
				{
					int r = pmap[y + ny][x + nx];
					surround[surround_index] = r;
					surround_index++;
				}
			}
		}

		part_temp_transition(this, i, x, y);
		part_pressure_transition(this, i, x, y);
		part_conduct_heat(this, i, x, y, surround);
	});
			//the main particle loop function, goes over all particles.

	for (int i = 0; i <= parts_lastActiveIndex; i++)
	{
		Particle* part = &parts[i];
		int t = part->type;
		if (t == 0)
			continue;
		Element* element = &elements[t];

		int x = (int)(part->x + 0.5f);
		int y = (int)(part->y + 0.5f);
		int cellX = x / CELL;
		int cellY = y / CELL;

		if(elementRecount)
			elementCount[t]++;

		//this kills any particle out of the screen, or in a wall where it isn't supposed to go
		if (bmap[y/CELL][x/CELL]==WL_DETECT && emap[y/CELL][x/CELL]<8)
			set_emap(x/CELL, y/CELL);

		float& cellVelocityX = velocityx(x, y);
		float& cellVelocityY = velocityy(x, y);
		float& cellPressure = pressure(x, y);

		//adding to velocity from the particle's velocity
		cellVelocityX = cellVelocityX * element->AirLoss + element->AirDrag * part->vx;
		cellVelocityY = cellVelocityY * element->AirLoss + element->AirDrag * part->vy;

		if (element->HotAir)
		{
			if (t==PT_GAS||t==PT_NBLE)
			{
				if (cellPressure < 3.5f)
					cellPressure += element->HotAir * (3.5f - cellPressure);
				if (y+CELL<YRES && pressure(x, y + CELL) < 3.5f)
					pressure(x, y + CELL) += element->HotAir * (3.5f - pressure(x, y + CELL));
				if (x+CELL<XRES)
				{
					if (pressure(x + CELL, y) < 3.5f)
						pressure(x + CELL, y) += element->HotAir * (3.5f - pressure(x + CELL, y));
					if (y+CELL<YRES && pressure(x + CELL, y + CELL) < 3.5f)
						pressure(x + CELL, y + CELL) += element->HotAir * (3.5f - pressure(x + CELL, y + CELL));
				}
			}
			else//add the hotair variable to the pressure map, like black hole, or white hole.
			{
				cellPressure += element->HotAir;
				if (y+CELL<YRES)
					pressure(x, y + CELL) += element->HotAir;
				if (x+CELL<XRES)
				{
					pressure(x + CELL, y) += element->HotAir;
					if (y+CELL<YRES)
						pressure(x + CELL, y + CELL) += element->HotAir;
				}
			}
		}

		float pGravX;
		float pGravY;

		if (element->Gravity || !(element->Properties & TYPE_SOLID))
		{
			get_gravity_at(this, x, y, element->Gravity, &pGravX, &pGravY);

			//Get some gravity from the gravity map
			if (t==PT_ANAR)
			{
				// perhaps we should have a ptypes variable for this
				pGravX -= gravx[(y/CELL)*(XRES/CELL)+(x/CELL)];
				pGravY -= gravy[(y/CELL)*(XRES/CELL)+(x/CELL)];
			}
			else if(t!=PT_STKM && t!=PT_STKM2 && t!=PT_FIGH && !(element->Properties & TYPE_SOLID))
			{
				pGravX += gravx[(y/CELL)*(XRES/CELL)+(x/CELL)];
				pGravY += gravy[(y/CELL)*(XRES/CELL)+(x/CELL)];
			}
		}
		else
			pGravX = pGravY = 0;
		//velocity updates for the particle
		if (t != PT_SPNG || !(part->flags&FLAG_MOVABLE))
		{
			part->vx *= element->Loss;
			part->vy *= element->Loss;
		}

		//particle gets velocity from the vx and vy maps
		part->vx += element->Advection * cellVelocityX + pGravX;
		part->vy += element->Advection * cellVelocityY + pGravY;

		if (element->Diffusion != 0)
		{
#ifdef REALISTIC
			part->vx += 0.05*sqrtf(part->temp)*element->Diffusion*(rand()/(0.5f*RAND_MAX)-1.0f);
			part->vy += 0.05*sqrtf(part->temp)*element->Diffusion*(rand()/(0.5f*RAND_MAX)-1.0f);
#else
			part->vx += element->Diffusion * (rand() / (0.5f * RAND_MAX) - 1.0f);
			part->vy += element->Diffusion * (rand() / (0.5f * RAND_MAX) - 1.0f);
#endif
		}

		int surround[8];
		int surround_index = 0;
		int surround_space = 0;
		int nt = 0;//if nt is greater than 1 after this, then there is a particle around the current particle, that is NOT the current particle's type, for water movement.
		for (int nx=-1; nx<2; nx++)
		{
			for (int ny=-1; ny<2; ny++)
			{
				if (nx||ny)
				{
					int r = pmap[y + ny][x + nx];
					surround[surround_index] = r;
					surround_index++;
					if (!(r&0xFF))
						surround_space++;//there is empty space
					if ((r&0xFF)!=t)
						nt++;//there is nothing or a different particle
				}
			}
		}

		part_conduct_heat(this, i, x, y, surround);

		if (t==PT_LIFE)
		{
			part->temp = restrict_flt(part->temp-50.0f, MIN_TEMP, MAX_TEMP);
			//ISGOL=1;//means there is a life particle on screen
		}
		if (t==PT_WIRE)
		{
			//wire_placed = 1;
		}
		//spark updates from walls
		if ((element->Properties & PROP_CONDUCTS) || t == PT_SPRK)
		{
			int nx = x % CELL;
			if (nx == 0)
				nx = x/CELL - 1;
			else if (nx == CELL-1)
				nx = x/CELL + 1;
			else
				nx = x/CELL;
			int ny = y % CELL;
			if (ny == 0)
				ny = y/CELL - 1;
			else if (ny == CELL-1)
				ny = y/CELL + 1;
			else
				ny = y/CELL;
			if (nx>=0 && ny>=0 && nx<XRES/CELL && ny<YRES/CELL)
			{
				if (t!=PT_SPRK)
				{
					if (emap[ny][nx]==12 && !part->life)
					{
						part_change_type(i,x,y,PT_SPRK);
						part->life = 4;
						part->ctype = t;
						t = PT_SPRK;
					}
				}
				else if (bmap[ny][nx]==WL_DETECT || bmap[ny][nx]==WL_EWALL || bmap[ny][nx]==WL_ALLOWLIQUID || bmap[ny][nx]==WL_WALLELEC || bmap[ny][nx]==WL_ALLOWALLELEC || bmap[ny][nx]==WL_EHOLE)
					set_emap(nx, ny);
			}
		}

		//the basic explosion, from the .explosive variable
		if ((element->Explosive&2) && cellPressure>2.5f)
		{
			part->life = rand()%80+180;
			part->temp = restrict_flt(elements[PT_FIRE].Temperature + (element->Flammable/2), MIN_TEMP, MAX_TEMP);
			t = PT_FIRE;
			part_change_type(i,x,y,t);
			cellPressure += 0.25f * CFDS;
		}

		//call the particle update function, if there is one
		if (element->Update && lua_el_mode[t] != 2)
		{
			if (element->Update(this, i, x, y, surround_space, nt, parts, pmap))
				continue;
			else if (t == PT_WARP)
			{
				// Warp does some movement in its update func, update variables to avoid incorrect data in pmap
				x = (int)(part->x+0.5f);
				y = (int)(part->y+0.5f);
			}
		}
		if(lua_el_mode[t])
		{
			if(luacon_elementReplacement(this, i, x, y, surround_space, nt, parts, pmap))
				continue;
			// Need to update variables, in case they've been changed by Lua
			x = (int)(part->x+0.5f);
			y = (int)(part->y+0.5f);
		}

		if (legacy_enable)//if heat sim is off
			Element::legacyUpdate(this, i,x,y,surround_space,nt, parts, pmap);

		if (part->type == PT_NONE)//if its dead, skip to next particle
			continue;

		part_movement(this, i, x, y, nt, surround_space);
	}
}

int part_pressure_transition(Simulation* sim, int i, int x, int y)
{
	Particle* part = &sim->parts[i];
	int t = part->type;
	Element* element = &sim->elements[t];

	float cellPressure = sim->pressure(x, y);
	bool typeChanged = true;
	float grav_pressure = 4 * fabs(sim->gravy[(y/CELL)*(XRES/CELL)+(x/CELL)])+fabs(sim->gravx[(y/CELL)*(XRES/CELL)+(x/CELL)]);
	float pressure = std::max<float>(grav_pressure * 4, cellPressure);

		// particle type change due to high pressure
	if (pressure > element->HighPressure && element->HighPressureTransition > -1)
	{
		if (element->HighPressureTransition != PT_NUM)
			t = element->HighPressureTransition;
		else if (t == PT_BMTL)
		{
			if (pressure > 2.5f || (pressure > 1.0f && part->tmp == 1))
				t = PT_BRMT;
			else
				typeChanged = false;
		}
		else
			typeChanged = false;
	}
		// particle type change due to low pressure
	else if (cellPressure < element->LowPressure && element->LowPressureTransition > -1)
	{
		if (element->LowPressureTransition != PT_NUM)
			t = element->LowPressureTransition;
		else
			typeChanged = false;
	}
	else
		typeChanged = false;
		// particle type change occurred
	if (typeChanged)
	{
		if (t == PT_NONE)
		{
			sim->kill_part(i);
			return t;
		}
		part->life = 0;
		sim->part_change_type(i,x,y,t);
		if (t == PT_FIRE)
			part->life = rand() % 50 + 120;
	}

	return t;
}

int part_temp_transition(Simulation* sim, int i, int x, int y)
{
	Particle* part = &sim->parts[i];
	int t = part->type;
	Element* element = &sim->elements[t];

	float temperature = part->temp;

	float ctemph = temperature;
	float ctempl = temperature;

	float cellPressure = sim->pressure(x, y);

	// change boiling point with pressure
	if ((element->State == ST_LIQUID &&
		element->HighTemperatureTransition > -1 &&element->HighTemperatureTransition < PT_NUM &&
		sim->elements[element->HighTemperatureTransition].State == ST_GAS) ||
		t == PT_LNTG || t == PT_SLTW)
	{
		ctemph -= 2.0f * cellPressure;
	}
	else if ((element->State == ST_GAS &&
		element->LowTemperatureTransition > -1 && element->LowTemperatureTransition < PT_NUM &&
		sim->elements[element->LowTemperatureTransition].State == ST_LIQUID) ||
		t == PT_WTRV)
	{
		ctempl -= 2.0f * cellPressure;
	}

	bool typeChanged = true;

	//A fix for ice with ctype = 0
	if ((t==PT_ICEI || t==PT_SNOW) && (part->ctype==0 || part->ctype>=PT_NUM || part->ctype==PT_ICEI || part->ctype==PT_SNOW))
		part->ctype = PT_WATR;

	// particle type change due to high temperature
	if (ctemph > element->HighTemperature && element->HighTemperatureTransition > -1)
	{
#ifdef REALISTIC
		float pt = part->temp;
		float dbt = ctempl - pt;
		if (element->HighTemperatureTransition!=PT_NUM)
		{
			if (platent[t] <= (c_heat - (element->HighTemperature - dbt)*c_Cm))
			{
				pt = (c_heat - platent[t])/c_Cm;
				t = element->HighTemperatureTransition;
			}
			else
			{
				part->temp = restrict_flt(element->HighTemperature - dbt, MIN_TEMP, MAX_TEMP);
				typeChanged = false;
			}
		}
#else
		if (element->HighTemperatureTransition != PT_NUM)
			t = element->HighTemperatureTransition;
#endif
		else if (t==PT_ICEI || t==PT_SNOW)
		{
			if (part->ctype < PT_NUM && part->ctype != t)
			{
				if (sim->elements[part->ctype].LowTemperatureTransition == t && temperature <= sim->elements[part->ctype].LowTemperature)
					typeChanged = false;
				else
				{
#ifdef REALISTIC
					//One ice table value for all it'typeChanged kinds
					if (platent[t] <= (c_heat - (sim->elements[part->ctype].LowTemperature - dbt)*c_Cm))
					{
						pt = (c_heat - platent[t])/c_Cm;
						t = part->ctype;
						part->ctype = PT_NONE;
						part->life = 0;
					}
					else
					{
						part->temp = restrict_flt(sim->elements[part->ctype].LowTemperature - dbt, MIN_TEMP, MAX_TEMP);
						typeChanged = false;
					}
#else
					t = part->ctype;
					part->ctype = PT_NONE;
					part->life = 0;
#endif
				}
			}
			else
				typeChanged = false;
		}
		else if (t==PT_SLTW) {
#ifdef REALISTIC
			if (platent[t] <= (c_heat - (element->HighTemperature - dbt)*c_Cm))
			{
				pt = (c_heat - platent[t])/c_Cm;

				if (rand()%4==0) t = PT_SALT;
				else t = PT_WTRV;
			}
			else
			{
				part->temp = restrict_flt(element->HighTemperature - dbt, MIN_TEMP, MAX_TEMP);
				typeChanged = false;
			}
#else
			if (rand()%4==0) t = PT_SALT;
			else t = PT_WTRV;
#endif
		}
		else typeChanged = false;
	} else if (ctempl < element->LowTemperature && element->LowTemperatureTransition > -1) {
		// particle type change due to low temperature
#ifdef REALISTIC
		float dbt = ctempl - temperature;
		if (element->LowTemperatureTransition!=PT_NUM)
		{
			if (platent[element->LowTemperatureTransition] >= (c_heat - (element->LowTemperature - dbt)*c_Cm))
			{
				temperature = (c_heat + platent[element->LowTemperatureTransition])/c_Cm;
				t = element->LowTemperatureTransition;
			}
			else
			{
				part->temp = restrict_flt(element->LowTemperature - dbt, MIN_TEMP, MAX_TEMP);
				typeChanged = false;
			}
		}
#else
		if (element->LowTemperatureTransition!=PT_NUM)
			t = element->LowTemperatureTransition;
#endif
		else if (t==PT_WTRV)
		{
			if (temperature < 273.0f)
				t = PT_RIME;
			else
				t = PT_DSTW;
		}
		else if (t==PT_LAVA) {
			if (part->ctype>0 && part->ctype<PT_NUM && part->ctype!=PT_LAVA) {
				if (part->ctype==PT_THRM&&temperature >= sim->elements[PT_BMTL].HighTemperature)
					typeChanged = false;
				else if ((part->ctype == PT_VIBR || part->ctype == PT_BVBR) && temperature >= 273.15f)
					typeChanged = false;
				else if (part->ctype==PT_TUGN && temperature > 3695.0f)
					typeChanged = false;
				else if (sim->elements[part->ctype].HighTemperatureTransition == PT_LAVA && temperature >= sim->elements[part->ctype].HighTemperature)
				{
					if (temperature>=sim->elements[part->ctype].HighTemperature)
						typeChanged = false;
				}
				else if (temperature >= 973.0f)
					typeChanged = false; // freezing point for lava with any other (not listed in ptransitions as turning into lava) ctype
				if (typeChanged)
				{
					t = part->ctype;
					part->ctype = PT_NONE;
					if (t==PT_THRM)
					{
						part->tmp = 0;
						t = PT_BMTL;
					}
					if (t==PT_PLUT)
					{
						part->tmp = 0;
						t = PT_LAVA;
					}
				}
			}
			else if (temperature < 973.0f)
				t = PT_STNE;
			else
				typeChanged = false;
		}
		else
			typeChanged = false;
	}
	else
		typeChanged = false;
#ifdef REALISTIC
	pt = restrict_flt(pt, MIN_TEMP, MAX_TEMP);
	for (j=0; j<8; j++)
	{
		parts[surround_hconduct[j]].temp = pt;
	}
#endif
	if (typeChanged)
	{ // particle type change occurred

		if (t == PT_NONE)
		{
			sim->kill_part(i);
			return t;
		}

		Element* newElement = &sim->elements[t];

		if (part->type != PT_ICEI && part->ctype != PT_FRZW)
			part->life = 0;

		if (t == PT_ICEI || t == PT_LAVA || t == PT_SNOW)
			part->ctype = part->type;
		if (newElement->State == ST_GAS && element->State != ST_GAS)
			sim->pressure(x, y) += 0.50f;
		if (t == PT_FIRE || t == PT_PLSM || t == PT_CFLM)
			part->life = rand() % 50 + 120;
		else if (t == PT_LAVA)
		{
			if (part->ctype == PT_BRMT)
				part->ctype = PT_BMTL;
			else if (part->ctype == PT_SAND)
				part->ctype = PT_GLAS;
			else if (part->ctype == PT_BGLA)
				part->ctype = PT_GLAS;
			else if (part->ctype == PT_PQRT)
				part->ctype = PT_QRTZ;
			part->life = rand() % 120 + 240;
		}
		sim->part_change_type(i,x,y,t);
	}

	temperature = part->temp;
	if (t==PT_LAVA) {
		part->life = restrict_flt((part->temp-700)/7, 0.0f, 400.0f);
		if (part->ctype==PT_THRM&&part->tmp>0)
		{
			part->tmp--;
			part->temp = 3500;
		}
		if (part->ctype==PT_PLUT&&part->tmp>0)
		{
			part->tmp--;
			part->temp = MAX_TEMP;
		}
	}

	return t;
}

float part_conduct_heat(Simulation* sim, int i, int x, int y, int surround[8])
{
	Particle* part = &sim->parts[i];
	int t = part->type;
	Element* element = &sim->elements[t];

	if (!sim->legacy_enable)
	{
		if (y >= 2 && y < YRES + 2 && (element->Properties&TYPE_LIQUID) && (t!=PT_GEL || (part->tmp * 2.55f) > (1 + rand() % 255))) //some heat convection for liquids
		{
			int r = sim->pmap[y-2][x];
			int otherType = r & 0xFF;
			int otherIndex = r >> 8;
			if (r != 0 && part->type == otherType && part->temp > sim->parts[otherIndex].temp)
				std::swap(part->temp, sim->parts[otherIndex].temp);
		}
	}

	float temperature = 0;
	bool conduct = true;
	if (t == PT_HSWC && part->life != 10)
		conduct = false;
	else if (t == PT_GEL)
	{
		float gel_scale = part->tmp * 2.55f;
		conduct = (element->HeatConduct * gel_scale) > (rand() % 250);
	}
	else
		conduct = element->HeatConduct > (rand() % 250);

			//heat transfer code
#ifdef REALISTIC
	if (t&&(t!=PT_HSWC||part->life==10)&&(element->HeatConduct*gel_scale))
#else
	if (conduct)
#endif
	{
		int h_count = 1;
		float c_Cm = 0.0f;
		float c_heat = 0.0f;
		if (sim->aheat_enable && !(element->Properties & PROP_NOAMBHEAT))
		{
#ifdef REALISTIC
			c_heat = part->temp*96.645/element->HeatConduct*gel_scale*fabs(element->Weight) + hv[y/CELL][x/CELL]*100*(cellPressure+273.15f)/256;
			c_Cm = 96.645/element->HeatConduct*gel_scale*fabs(element->Weight)  + 100*(cellPressure+273.15f)/256;
			pt = c_heat/c_Cm;
			pt = restrict_flt(pt, -MAX_TEMP+MIN_TEMP, MAX_TEMP-MIN_TEMP);
			part->temp = pt;
			//Pressure increase from heat (temporary)
			cellPressure += (pt-hv[y/CELL][x/CELL])*0.004;
			hv[y/CELL][x/CELL] = pt;
#else
			float& aheat = sim->ambientheat(x, y);
			c_heat = (aheat - part->temp)*0.04;
			c_heat = restrict_flt(c_heat, -MAX_TEMP+MIN_TEMP, MAX_TEMP-MIN_TEMP);
			part->temp += c_heat;
			aheat -= c_heat;
#endif
		}
		c_heat = part->temp;
		c_Cm = 0.0f;

		int surround_hconduct[8];
		for (int j = 0; j < 8; j++)
		{
			surround_hconduct[j] = i;
			int r = surround[j];
			if (!r)
				continue;
			int rt = r&0xFF;
			if (rt && sim->elements[rt].HeatConduct && (rt!=PT_HSWC || sim->parts[r>>8].life==10)
				&& (t!=PT_FILT||(rt!=PT_BRAY&&rt!=PT_BIZR&&rt!=PT_BIZRG))
				&& (rt!=PT_FILT||(t!=PT_BRAY&&t!=PT_PHOT&&t!=PT_BIZR&&t!=PT_BIZRG))
				&& (t!=PT_ELEC||rt!=PT_DEUT)
				&& (t!=PT_DEUT||rt!=PT_ELEC))
			{
				surround_hconduct[j] = r>>8;
#ifdef REALISTIC
				if (rt==PT_GEL)
					gel_scale = parts[r>>8].tmp*2.55f;
				else gel_scale = 1.0f;

				c_heat += parts[r>>8].temp*96.645/elements[rt].HeatConduct*gel_scale*fabs(elements[rt].Weight);
				c_Cm += 96.645/elements[rt].HeatConduct*gel_scale*fabs(elements[rt].Weight);
#else
				c_heat += sim->parts[r>>8].temp;
#endif
				h_count++;
			}
		}
#ifdef REALISTIC
		if (t==PT_GEL)
			gel_scale = part->tmp*2.55f;
		else gel_scale = 1.0f;

		if (t == PT_PHOT)
			pt = (c_heat+part->temp*96.645)/(c_Cm+96.645);
		else
			pt = (c_heat+part->temp*96.645/element->HeatConduct*gel_scale*fabs(element->Weight))/(c_Cm+96.645/element->HeatConduct*gel_scale*fabs(element->Weight));

		c_heat += part->temp*96.645/element->HeatConduct*gel_scale*fabs(element->Weight);
		c_Cm += 96.645/element->HeatConduct*gel_scale*fabs(element->Weight);
		part->temp = restrict_flt(pt, MIN_TEMP, MAX_TEMP);
#else
		temperature = c_heat / h_count;
		temperature = part->temp = restrict_flt(temperature, MIN_TEMP, MAX_TEMP);
		for (int j=0; j<8; j++)
		{
			sim->parts[surround_hconduct[j]].temp = temperature;
		}
#endif
	}
	else
		part->temp = restrict_flt(part->temp, MIN_TEMP, MAX_TEMP);

	return part->temp;
}

void part_movement(Simulation* sim, int i, int x, int y, int nt, int surround_space)
{
	Particle* part = &sim->parts[i];
	int t = part->type;
	Element* element = &sim->elements[t];

	if (!part->vx && !part->vy)//if its not moving, skip to next particle, movement code is next
		return;

	float mv = std::max<float>(fabsf(part->vx), fabsf(part->vy));
	int clear_x;
	int clear_y;
	int fin_x;
	int fin_y;
	float clear_xf;
	float clear_yf;
	float fin_xf;
	float fin_yf;
	if (mv < ISTP)
	{
		clear_x = x;
		clear_y = y;
		clear_xf = part->x;
		clear_yf = part->y;
		fin_xf = clear_xf + part->vx;
		fin_yf = clear_yf + part->vy;
		fin_x = (int)(fin_xf+0.5f);
		fin_y = (int)(fin_yf+0.5f);
	}
	else
	{
		// interpolate to see if there is anything in the way
		float dx = part->vx*ISTP/mv;
		float dy = part->vy*ISTP/mv;
		fin_xf = part->x;
		fin_yf = part->y;
		while (1)
		{
			mv -= ISTP;
			fin_xf += dx;
			fin_yf += dy;
			fin_x = (int)(fin_xf+0.5f);
			fin_y = (int)(fin_yf+0.5f);
			if (mv <= 0.0f)
			{
				// nothing found
				fin_xf = part->x + part->vx;
				fin_yf = part->y + part->vy;
				fin_x = (int)(fin_xf+0.5f);
				fin_y = (int)(fin_yf+0.5f);
				clear_xf = fin_xf-dx;
				clear_yf = fin_yf-dy;
				clear_x = (int)(clear_xf+0.5f);
				clear_y = (int)(clear_yf+0.5f);
				break;
			}
			if (fin_x<CELL || fin_y<CELL || fin_x>=XRES-CELL || fin_y>=YRES-CELL || sim->pmap[fin_y][fin_x] || (sim->bmap[fin_y/CELL][fin_x/CELL] && (sim->bmap[fin_y/CELL][fin_x/CELL]==WL_DESTROYALL || !sim->eval_move(t,fin_x,fin_y,NULL))))
			{
				// found an obstacle
				clear_xf = fin_xf-dx;
				clear_yf = fin_yf-dy;
				clear_x = (int)(clear_xf+0.5f);
				clear_y = (int)(clear_yf+0.5f);
				break;
			}
			if (sim->bmap[fin_y/CELL][fin_x/CELL]==WL_DETECT && sim->emap[fin_y/CELL][fin_x/CELL]<8)
				sim->set_emap(fin_x/CELL, fin_y/CELL);
		}
	}

	bool stagnant = (part->flags & FLAG_STAGNANT) != 0;
	part->flags &= ~FLAG_STAGNANT;

	if (t==PT_STKM || t==PT_STKM2 || t==PT_FIGH)
	{
		int nx, ny;
		//head movement, let head pass through anything
		part->x += part->vx;
		part->y += part->vy;
		nx = (int)((float)part->x+0.5f);
		ny = (int)((float)part->y+0.5f);
		if (ny!=y || nx!=x)
		{
			if ((sim->pmap[y][x]>>8)==i) sim->pmap[y][x] = 0;
			else if ((sim->photons[y][x]>>8)==i) sim->photons[y][x] = 0;
			if (nx<CELL || nx>=XRES-CELL || ny<CELL || ny>=YRES-CELL)
			{
				sim->kill_part(i);
				return;
			}
			if (element->Properties & TYPE_ENERGY)
				sim->photons[ny][nx] = t|(i<<8);
			else if (t)
				sim->pmap[ny][nx] = t|(i<<8);
		}
	}
	else if (element->Properties & TYPE_ENERGY)
	{
		if (t == PT_PHOT) {
			if (part->flags&FLAG_SKIPMOVE)
			{
				part->flags &= ~FLAG_SKIPMOVE;
				return;
			}

			int rt = sim->pmap[fin_y][fin_x] & 0xFF;
			int lt = sim->pmap[y][x] & 0xFF;

			int r = sim->eval_move(PT_PHOT, fin_x, fin_y, NULL);
			if (((rt==PT_GLAS && lt!=PT_GLAS) || (rt!=PT_GLAS && lt==PT_GLAS)) && r) {
				float nrx;
				float nry;

				if (!sim->get_normal_interp(REFRACT|t, part->x, part->y, part->vx, part->vy, &nrx, &nry)) {
					sim->kill_part(i);
					return;
				}

				r = sim->get_wavelength_bin(&part->ctype);
				if (r == -1) {
					sim->kill_part(i);
					return;
				}
				int nn = GLASS_IOR - GLASS_DISP*(r-15)/15.0f;
				nn *= nn;
				nrx = -nrx;
				nry = -nry;
				if (rt==PT_GLAS && lt!=PT_GLAS)
					nn = 1.0f/nn;
				int ct1 = part->vx*nrx + part->vy*nry;
				int ct2 = 1.0f - (nn*nn)*(1.0f-(ct1*ct1));
				if (ct2 < 0.0f) {
					// total internal reflection
					part->vx -= 2.0f*ct1*nrx;
					part->vy -= 2.0f*ct1*nry;
					fin_xf = part->x;
					fin_yf = part->y;
					fin_x = x;
					fin_y = y;
				} else {
					// refraction
					ct2 = sqrtf(ct2);
					ct2 = ct2 - nn*ct1;
					part->vx = nn*part->vx + ct2*nrx;
					part->vy = nn*part->vy + ct2*nry;
				}
			}
		}
		if (stagnant)//FLAG_STAGNANT set, was reflected on previous frame
		{
			// cast coords as int then back to float for compatibility with existing saves
			if (!sim->do_move(i, x, y, (float)fin_x, (float)fin_y) && part->type) {
				sim->kill_part(i);
				return;
			}
		}
		else if (!sim->do_move(i, x, y, fin_xf, fin_yf))
		{
			if (part->type == PT_NONE)
				return;
			// reflection
			part->flags |= FLAG_STAGNANT;
			if (t==PT_NEUT && 100>(rand()%1000))
			{
				sim->kill_part(i);
				return;
			}
			int r = sim->pmap[fin_y][fin_x];

			if (((r&0xFF)==PT_PIPE || (r&0xFF) == PT_PPIP) && !(sim->parts[r>>8].tmp&0xFF))
			{
				sim->parts[r>>8].tmp =  (sim->parts[r>>8].tmp&~0xFF) | part->type;
				sim->parts[r>>8].temp = part->temp;
				sim->parts[r>>8].tmp2 = part->life;
				sim->parts[r>>8].pavg[0] = part->tmp;
				sim->parts[r>>8].pavg[1] = part->ctype;
				sim->kill_part(i);
				return;
			}

			// this should be replaced with a particle type attribute ("photwl" or something)
			if ((r & 0xFF) == PT_PSCN) part->ctype  = 0x00000000;
			if ((r & 0xFF) == PT_NSCN) part->ctype  = 0x00000000;
			if ((r & 0xFF) == PT_SPRK) part->ctype  = 0x00000000;
			if ((r & 0xFF) == PT_COAL) part->ctype  = 0x00000000;
			if ((r & 0xFF) == PT_BCOL) part->ctype  = 0x00000000;
			if ((r & 0xFF) == PT_PLEX) part->ctype &= 0x1F00003E;
			if ((r & 0xFF) == PT_NITR) part->ctype &= 0x0007C000;
			if ((r & 0xFF) == PT_NBLE) part->ctype &= 0x3FFF8000;
			if ((r & 0xFF) == PT_LAVA) part->ctype &= 0x3FF00000;
			if ((r & 0xFF) == PT_ACID) part->ctype &= 0x1FE001FE;
			if ((r & 0xFF) == PT_DUST) part->ctype &= 0x3FFFFFC0;
			if ((r & 0xFF) == PT_SNOW) part->ctype &= 0x03FFFFFF;
			if ((r & 0xFF) == PT_GOO)  part->ctype &= 0x3FFAAA00;
			if ((r & 0xFF) == PT_PLNT) part->ctype &= 0x0007C000;
			if ((r & 0xFF) == PT_PLUT) part->ctype &= 0x001FCE00;
			if ((r & 0xFF) == PT_URAN) part->ctype &= 0x003FC000;

			float nrx;
			float nry;
			if (sim->get_normal_interp(t, part->x, part->y, part->vx, part->vy, &nrx, &nry)) {
				float dp = nrx*part->vx + nry*part->vy;
				part->vx -= 2.0f*dp*nrx;
				part->vy -= 2.0f*dp*nry;
				// leave the actual movement until next frame so that reflection of fast particles and refraction happen correctly
			} else {
				if (t!=PT_NEUT)
					sim->kill_part(i);
				return;
			}
			if (!(part->ctype&0x3FFFFFFF) && t == PT_PHOT) {
				sim->kill_part(i);
				return;
			}
		}
	}
	else if (element->Falldown==0)
	{
		// gasses and solids (but not powders)
		if (!sim->do_move(i, x, y, fin_xf, fin_yf))
		{
			if (part->type == PT_NONE)
				return;
			// can't move there, so bounce off
			// TODO
			// TODO: Work out what previous TODO was for
			if (fin_x>x+ISTP) fin_x=x+ISTP;
			if (fin_x<x-ISTP) fin_x=x-ISTP;
			if (fin_y>y+ISTP) fin_y=y+ISTP;
			if (fin_y<y-ISTP) fin_y=y-ISTP;
			if (sim->do_move(i, x, y, 0.25f+(float)(2*x-fin_x), 0.25f+fin_y))
			{
				part->vx *= element->Collision;
			}
			else if (sim->do_move(i, x, y, 0.25f+fin_x, 0.25f+(float)(2*y-fin_y)))
			{
				part->vy *= element->Collision;
			}
			else
			{
				part->vx *= element->Collision;
				part->vy *= element->Collision;
			}
		}
	}
	else
	{
		if (sim->water_equal_test && element->Falldown == 2 && 1>= rand()%400)//checking stagnant is cool, but then it doesn't update when you change it later.
		{
			if (!sim->flood_water(x,y,i,y, part->flags&FLAG_WATEREQUAL))
				return;
		}

		// liquids and powders
		if (!sim->do_move(i, x, y, fin_xf, fin_yf))
		{
			if (part->type == PT_NONE)
				return;
			if (fin_x!=x && sim->do_move(i, x, y, fin_xf, clear_yf))
			{
				part->vx *= element->Collision;
				part->vy *= element->Collision;
			}
			else if (fin_y!=y && sim->do_move(i, x, y, clear_xf, fin_yf))
			{
				part->vx *= element->Collision;
				part->vy *= element->Collision;
			}
			else
			{
				float pGravX;
				float pGravY;

				get_gravity_at(sim, x, y, element->Gravity, &pGravX, &pGravY);

				int r = (rand()%2)*2-1;
				// allow diagonal movement if target position is blocked
				// but no point trying this if particle is stuck in a block of identical particles
				if ((clear_x!=x || clear_y!=y || nt || surround_space) &&
					(fabsf(part->vx)>0.01f || fabsf(part->vy)>0.01f))
				{
					float dx = part->vx - part->vy * r;
					float dy = part->vy + part->vx * r;

					mv = std::max<float>(fabsf(dx), fabsf(dy));

					dx /= mv;
					dy /= mv;
					// Try to move to one of the diagonals
					if (sim->do_move(i, x, y, clear_xf+dx, clear_yf+dy))
					{
						part->vx *= element->Collision;
						part->vy *= element->Collision;
						return;
					}
					std::swap(dx, dy);
					dx *= r;
					dy *= -r;
					// Now try the other one
					if (sim->do_move(i, x, y, clear_xf+dx, clear_yf+dy))
					{
						part->vx *= element->Collision;
						part->vy *= element->Collision;
						return;
					}
				}
				if (element->Falldown>1 && !sim->grav->ngrav_enable && sim->gravityMode==0 && part->vy>fabsf(part->vx))
				{
					int move_result = 0;

					int rt;
					// stagnant is true if FLAG_STAGNANT was set for this particle in previous frame
					if (!stagnant || nt) //nt is if there is an something else besides the current particle type, around the particle
						rt = 30;//slight less water lag, although it changes how it moves a lot
					else
						rt = 10;

					if (t==PT_GEL)
						rt = part->tmp*0.20f+5.0f;

					int nx;
					int ny;

					for (int j=clear_x+r; j>=0 && j>=clear_x-rt && j<clear_x+rt && j<XRES; j+=r)
					{
						if (((sim->pmap[fin_y][j]&0xFF)!=t || sim->bmap[fin_y/CELL][j/CELL])
							&& (move_result = sim->do_move(i, x, y, (float)j, fin_yf)))
						{
							nx = (int)(part->x+0.5f);
							ny = (int)(part->y+0.5f);
							break;
						}
						if (fin_y!=clear_y && ((sim->pmap[clear_y][j]&0xFF)!=t || sim->bmap[clear_y/CELL][j/CELL])
							&& (move_result = sim->do_move(i, x, y, (float)j, clear_yf)))
						{
							nx = (int)(part->x+0.5f);
							ny = (int)(part->y+0.5f);
							break;
						}
						if ((sim->pmap[clear_y][j]&0xFF)!=t || (sim->bmap[clear_y/CELL][j/CELL] && sim->bmap[clear_y/CELL][j/CELL]!=WL_STREAM))
							break;
					}
					if (part->vy>0)
						r = 1;
					else
						r = -1;
					if (move_result == 1)
					{
						for (int j=ny+r; j>=0 && j<YRES && j>=ny-rt && j<ny+rt; j+=r)
						{
							if (((sim->pmap[j][nx]&0xFF)!=t || sim->bmap[j/CELL][nx/CELL]) && sim->do_move(i, nx, ny, (float)nx, (float)j))
								break;
							if ((sim->pmap[j][nx]&255)!=t || (sim->bmap[j/CELL][nx/CELL] && sim->bmap[j/CELL][nx/CELL]!=WL_STREAM))
								break;
						}
					}
					else if (move_result == -1) {} // particle is out of bounds
					else if ((clear_x!=x||clear_y!=y) && sim->do_move(i, x, y, clear_xf, clear_yf)) {}
					else part->flags |= FLAG_STAGNANT;
					part->vx *= element->Collision;
					part->vy *= element->Collision;
				}
				else if (element->Falldown>1 && fabsf(pGravX*part->vx+pGravY*part->vy)>fabsf(pGravY*part->vx-pGravX*part->vy))
				{
					int move_result = 0;

					int rt;
					// stagnant is true if FLAG_STAGNANT was set for this particle in previous frame
					if (!stagnant || nt) //nt is if there is an something else besides the current particle type, around the particle
						rt = 30;//slight less water lag, although it changes how it moves a lot
					else
						rt = 10;
					float nxf = clear_xf;
					float nyf = clear_yf;
					int nx = (int) (nxf + 0.5f);
					int ny = (int) (nyf + 0.5f);

					float pGravX;
					float pGravY;
					float prev_pGravX;
					float prev_pGravY;
					float ptGrav = element->Gravity;

					for (int j=0;j<rt;j++)
					{
						get_gravity_at(sim, nx, ny, element->Gravity, &pGravX, &pGravY);

						pGravX += sim->gravx[(ny/CELL)*(XRES/CELL)+(nx/CELL)];
						pGravY += sim->gravy[(ny/CELL)*(XRES/CELL)+(nx/CELL)];
						if (fabsf(pGravY)>fabsf(pGravX))
							mv = fabsf(pGravY);
						else
							mv = fabsf(pGravX);
						if (mv<0.0001f) break;
						pGravX /= mv;
						pGravY /= mv;
						if (j)
						{
							nxf += r*(pGravY*2.0f-prev_pGravY);
							nyf += -r*(pGravX*2.0f-prev_pGravX);
						}
						else
						{
							nxf += r*pGravY;
							nyf += -r*pGravX;
						}
						prev_pGravX = pGravX;
						prev_pGravY = pGravY;
						nx = (int)(nxf+0.5f);
						ny = (int)(nyf+0.5f);
						if (nx<0 || ny<0 || nx>=XRES || ny >=YRES)
							break;
						if ((sim->pmap[ny][nx]&0xFF)!=t || sim->bmap[ny/CELL][nx/CELL])
						{
							move_result = sim->do_move(i, x, y, nxf, nyf);
							if (move_result != 0)
							{
								nx = (int)(part->x+0.5f);
								ny = (int)(part->y+0.5f);
								break;
							}
							if (sim->bmap[ny/CELL][nx/CELL]!=WL_STREAM)
								break;
						}
					}
					if (move_result == 1)
					{
						clear_x = nx;
						clear_y = ny;
						for (int j = 0; j < rt; j++)
						{
							get_gravity_at(sim, nx, ny, element->Gravity, &pGravX, &pGravY);

							pGravX += sim->gravx[(ny/CELL)*(XRES/CELL)+(nx/CELL)];
							pGravY += sim->gravy[(ny/CELL)*(XRES/CELL)+(nx/CELL)];
							if (fabsf(pGravY)>fabsf(pGravX))
								mv = fabsf(pGravY);
							else
								mv = fabsf(pGravX);
							if (mv<0.0001f) break;
							pGravX /= mv;
							pGravY /= mv;
							nxf += pGravX;
							nyf += pGravY;
							nx = (int)(nxf+0.5f);
							ny = (int)(nyf+0.5f);
							if (nx<0 || ny<0 || nx>=XRES || ny>=YRES)
								break;
							if ((sim->pmap[ny][nx]&0xFF)!=t || sim->bmap[ny/CELL][nx/CELL])
							{
								move_result = sim->do_move(i, clear_x, clear_y, nxf, nyf);
								if (move_result != 0 || sim->bmap[ny/CELL][nx/CELL]!=WL_STREAM)
									break;
							}
						}
					}
					else if (move_result == -1) {} // particle is out of bounds
					else if ((clear_x!=x||clear_y!=y) && sim->do_move(i, x, y, clear_xf, clear_yf)) {}
					else part->flags |= FLAG_STAGNANT;
					part->vx *= element->Collision;
					part->vy *= element->Collision;
				}
				else
				{
					// if interpolation was done, try moving to last clear position
					if ((clear_x!=x||clear_y!=y) && sim->do_move(i, x, y, clear_xf, clear_yf)) {}
					else part->flags |= FLAG_STAGNANT;
					part->vx *= element->Collision;
					part->vy *= element->Collision;
				}
			}
		}
	}
}

bool get_gravity_at(Simulation* sim, int x, int y, float elementGravity, float* gravX, float* gravY)
{
	//Gravity mode by Moach
	switch (sim->gravityMode)
	{
	default:
	case 0:
		*gravX = 0.0f;
		*gravY = elementGravity;
		break;
	case 1:
		*gravX = *gravY = 0.0f;
		break;
	case 2:
		{
			float pGravD = 0.01f - hypotf((x - XCNTR), (y - YCNTR));
			*gravX = elementGravity * ((float)(x - XCNTR) / pGravD);
			*gravY = elementGravity * ((float)(y - YCNTR) / pGravD);
		}
		break;
	}

	return true;
}