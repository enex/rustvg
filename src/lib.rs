use std::default::Default;

#[derive(Clone, Debug)]
pub struct Paint;

#[derive(Debug, Clone)]
pub enum Commands{
    MOVETO,
    LINETO,
    BEZIERTO,
    CLOSE,
    WINDING,
}

#[derive(Debug, Clone)]
pub enum ExpandFeatures{
    FILL =   0x01,
    STROKE = 0x02,
    CAPS =   0x04,
}

#[derive(Debug, Clone)]
pub struct Scissor {
	xform: [f32; 6],
	extent: [f32; 2],
}

#[derive(Clone, Debug)]
pub struct State{
    fill: Paint,
    stroke: Paint,
    stroke_width: f32,
    miter_limit: f32,
    line_join: f32,
    line_cap: f32,
    xform: [f32; 6],
    scissor: Scissor,
    font_size: f32,
    letter_spacing: f32,
    line_height: f32,
    font_blur: f32,
    text_align: isize,
    font_id: isize,
    alpha: f32,
}

impl Default for State{
    fn default() -> State{
        State{
            fill: Paint,//TODO: do correct
            stroke: Paint,
            stroke_width: 1.0,
            miter_limit: 10.0,
            line_cap: 0.0,
            line_join: 0.0,
            alpha: 1.0,
            xform: [1.,0.,0.,1.,0.,0.],

            scissor: Scissor{
                xform: [1.,0.,0.,1.,0.,0.],
                extent: [-1., -1.],
            },

            font_size: 16.0,
            letter_spacing: 0.0,
            line_height: 1.0,
            font_blur: 0.0,
            text_align: 0,
            font_id: 0,
        }
    }
}

#[derive(Debug, Default, Clone)]
pub struct Point{
    x: f32,
    y: f32,
    dx: f32,
    dy: f32,
    len: f32,
    dmx: f32,
    dmy: f32,
    flags: usize,
}

pub trait Params{
    fn render_create(&mut self);
}

/*
struct NVGparams {
	void* userPtr;
	int edgeAntiAlias;
	int (*renderCreate)(void* uptr);
	int (*renderCreateTexture)(void* uptr, int type, int w, int h, int imageFlags, const unsigned char* data);
	int (*renderDeleteTexture)(void* uptr, int image);
	int (*renderUpdateTexture)(void* uptr, int image, int x, int y, int w, int h, const unsigned char* data);
	int (*renderGetTextureSize)(void* uptr, int image, int* w, int* h);
	void (*renderViewport)(void* uptr, int width, int height);
	void (*renderCancel)(void* uptr);
	void (*renderFlush)(void* uptr);
	void (*renderFill)(void* uptr, NVGpaint* paint, NVGscissor* scissor, float fringe, const float* bounds, const NVGpath* paths, int npaths);
	void (*renderStroke)(void* uptr, NVGpaint* paint, NVGscissor* scissor, float fringe, float strokeWidth, const NVGpath* paths, int npaths);
	void (*renderTriangles)(void* uptr, NVGpaint* paint, NVGscissor* scissor, const NVGvertex* verts, int nverts);
	void (*renderDelete)(void* uptr);
};*/

#[derive(Debug, Default, Clone)]
pub struct Vertex{
    x: f32,
    y: f32,
    u: f32,
    v: f32,
}

#[derive(Debug, Default, Clone)]
pub struct Path{
    first: usize,
    count: usize,
    closed: bool,
    nbevel: usize,
    fill: Vertex,
    stroke: Vertex,
    winding: usize,
    convex: usize,
}

#[derive(Debug, Clone)]
pub struct PathCache{
    points: Vec<Point>,
    paths: Vec<Path>,
    verts: Vec<Vertex>,
    bounds: [f32; 4],
}

#[derive(Debug, Clone)]
pub struct Context{
    //params: Params,
    commands: Vec<f32>,
    command_x: f32,
    command_y: f32,
    states: Vec<State>,
    cache: PathCache,
    tess_tol: f32,
    dist_tol: f32,
    fringe_width: f32,
    device_px_ratio: f32,
    //fs: FONScontext,
    font_image: isize,
    alpha_blend: isize,
    draw_call_count: isize,
    fill_tri_count: isize,
    sroke_tri_count: isize,
    text_tri_count: isize,
}

impl Context{
    pub fn begin_frame(&mut self, width: isize, height: isize, device_pixel_ratio: f32, alpha_blend: isize){
        self.states.truncate(0);
        self.draw_call_count = 0;
        self.fill_tri_count = 0;
        self.sroke_tri_count = 0;
        self.text_tri_count = 0;
    }
    pub fn end_frame(&mut self){

    }
    /// save the current state
    pub fn save(&mut self){
        let s = self.states[self.states.len()-1].clone();
        self.states.push(s);
    }
    pub fn restore(&mut self){
        self.states.pop();
    }
    /// reset current state do default state
    pub fn reset(&mut self){
        unimplemented!()
    }

    // some setter and getter functions

    /// add a point nvg__addPoint
    fn add_point(&mut self, x: f32, y: f32, tp: usize){

    }

    fn teslate_bezier(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x3: f32, y3: f32, x4: f32, y4: f32, level: usize, tp: usize){
        if level > 10{
            return;
        }

        let x12 = (x1+x2)*0.5;
    	let y12 = (y1+y2)*0.5;
    	let x23 = (x2+x3)*0.5;
    	let y23 = (y2+y3)*0.5;
    	let x34 = (x3+x4)*0.5;
    	let y34 = (y3+y4)*0.5;
    	let x123 = (x12+x23)*0.5;
    	let y123 = (y12+y23)*0.5;

        let dx = x3 - x1;
    	let dy = y3 - y1;
    	let d2 = (((x2 - x4) * dy - (y2 - y4) * dx)).abs();
    	let d3 = (((x3 - x4) * dy - (y3 - y4) * dx)).abs();

        if ((d2 + d3)*(d2 + d3) < self.tess_tol * (dx*dx + dy*dy)) {
    		self.add_point(x4, y4, tp);
    		return;
    	}

        let x234 = (x23+x34)*0.5;
    	let y234 = (y23+y34)*0.5;
    	let x1234 = (x123+x234)*0.5;
    	let y1234 = (y123+y234)*0.5;

    	self.teslate_bezier(x1,y1, x12,y12, x123,y123, x1234,y1234, level+1, 0);
    	self.teslate_bezier(x1234,y1234, x234,y234, x34,y34, x4,y4, level+1, tp);
    }

    fn flatten_path(&mut self){
        let ref cache = self.cache;
            /*
        	struct NVGpoint* last;
        	struct NVGpoint* p0;
        	struct NVGpoint* p1;
        	struct NVGpoint* pts;
        	struct NVGpath* path;
        	int i, j;
        	float* cp1;
        	float* cp2;
        	float* p;
        	float area;

        	if (cache->npaths > 0)
        		return;

        	// Flatten
        	i = 0;
        	while (i < ctx->ncommands) {
        		int cmd = (int)ctx->commands[i];
        		switch (cmd) {
        		case NVG_MOVETO:
        			nvg__addPath(ctx);
        			p = &ctx->commands[i+1];
        			nvg__addPoint(ctx, p[0], p[1], NVG_PT_CORNER);
        			i += 3;
        			break;
        		case NVG_LINETO:
        			p = &ctx->commands[i+1];
        			nvg__addPoint(ctx, p[0], p[1], NVG_PT_CORNER);
        			i += 3;
        			break;
        		case NVG_BEZIERTO:
        			last = nvg__lastPoint(ctx);
        			if (last != NULL) {
        				cp1 = &ctx->commands[i+1];
        				cp2 = &ctx->commands[i+3];
        				p = &ctx->commands[i+5];
        				nvg__tesselateBezier(ctx, last->x,last->y, cp1[0],cp1[1], cp2[0],cp2[1], p[0],p[1], 0, NVG_PT_CORNER);
        			}
        			i += 7;
        			break;
        		case NVG_CLOSE:
        			nvg__closePath(ctx);
        			i++;
        			break;
        		case NVG_WINDING:
        			nvg__pathWinding(ctx, (int)ctx->commands[i+1]);
        			i += 2;
        			break;
        		default:
        			i++;
        		}
        	}

        	cache->bounds[0] = cache->bounds[1] = 1e6f;
        	cache->bounds[2] = cache->bounds[3] = -1e6f;

        	// Calculate the direction and length of line segments.
        	for (j = 0; j < cache->npaths; j++) {
        		path = &cache->paths[j];
        		pts = &cache->points[path->first];

        		// If the first and last points are the same, remove the last, mark as closed path.
        		p0 = &pts[path->count-1];
        		p1 = &pts[0];
        		if (nvg__ptEquals(p0->x,p0->y, p1->x,p1->y, ctx->distTol)) {
        			path->count--;
        			p0 = &pts[path->count-1];
        			path->closed = 1;
        		}

        		// Enforce winding.
        		if (path->count > 2) {
        			area = nvg__polyArea(pts, path->count);
        			if (path->winding == NVG_CCW && area < 0.0f)
        				nvg__polyReverse(pts, path->count);
        			if (path->winding == NVG_CW && area > 0.0f)
        				nvg__polyReverse(pts, path->count);
        		}

        		for(i = 0; i < path->count; ++i) {
        			// Calculate segment direction and length
        			p0->dx = p1->x - p0->x;
        			p0->dy = p1->y - p0->y;
        			p0->len = nvg__normalize(&p0->dx, &p0->dy);
        			// Update bounds
        			cache->bounds[0] = nvg__minf(cache->bounds[0], p0->x);
        			cache->bounds[1] = nvg__minf(cache->bounds[1], p0->y);
        			cache->bounds[2] = nvg__maxf(cache->bounds[2], p0->x);
        			cache->bounds[3] = nvg__maxf(cache->bounds[3], p0->y);
        			// Advance
        			p0 = p1++;
        		}
        	}
        }*/
    }
}

/// https://github.com/DISTRHO/DPF/blob/master/dgl/src/nanovg/nanovg.c
#[test]
fn it_works() {
}
