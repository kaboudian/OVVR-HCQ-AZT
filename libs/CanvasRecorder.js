function readOption( opt, def ){
    return (opt) ? opt : def ;
}

class CanvasStreamRecorder{
    constructor( canvas, opt={} ){
        this._canvas  = canvas ;
        this._stream = canvas.captureStream() ;
        this._recorder = new MediaRecorder( this.stream, 
            { mimType : 'video/mp4;codecs=vp9' }  ) ;
        this._chunks = [] ;
        this._filename = readOption(opt.fileName , this.canvas.id+".webm" ) ;
        this._filename = readOption(opt.filename , this._filename ) ;
        this._recording = false ;
        this.recorder.canvasRecorder  = this ;
        this.recorder.ondataavailable = this.handleData ;
        this._resetting = false ;
    }
    get canvas(){
        return this._canvas ;
    }

    get stream(){
        return this._stream ;
    }

    get options(){
        return this._options ; 
    } 

    get recorder(){
        return this._recorder ;
    }

    get chunks(){
        return this._chunks ;
    }

    get filename(){
        return this._filename ;
    }

    set filename(nf){
        this._filename = nf ;
    }

    get active(){
        return ( this.recorder.state != 'inactive' ) ;
    }
    get recording(){
        return this._recording ;
    }

    set recording(nv){
        if (nv){
            if (!this.recording){
                this.start() ;
            }
        }else{
            if (this.recording){
                this._recording = false ;
                this.stop() ;
            }
        }
    }
    start(){
        this.paused = false ;
        if ( this.active ){
            console.log('Already running!') ;
            return ;
        }
        console.log('Recording started!') ;
        return this.recorder.start() ;
    }

    stop(){
        if ( this.active ){
            console.log('Recording stoped!') ;
            return this.recorder.stop() ;
        }else
            return false ;
    }

    pause(){
        this.paused = true ;
        console.log('Recording paused!') ;
        return (this.recording) ? this.recorder.pause() : null ;
    }

    resume(){
        this.paused = false ;
        console.log('Recording resumed!') ;
        return (this.active) ? this.recorder.resume() : this.start() ;
    }

    reset(){
        this._chunks = [] ;
        this._resetting = true ;
        this.stop() ;
    }

    handleData(event){
        //console.log(this.canvasRecorder) ;
        if (event.data.size > 0) {
            this.canvasRecorder.chunks.push(event.data);
            this.canvasRecorder.download();
        } else {
            console.log('Nothing to write to disk!') ;
        }
    }

    download(){
        //if ( this._resetting ){
        //    if (this.recording){
        //        this.start() ;
        //    }
        //    this._resetting = false ;
        //    return ;
        //}
        var blob = new Blob(this.chunks, {
            type: "video/mp4"
        });

        var url = URL.createObjectURL(blob);
        var a = document.createElement("a");
        a.style = "display: none";
        a.href = url;
        a.download = this.filename ;
        document.body.appendChild(a);
        a.click();
        window.URL.revokeObjectURL(url);
        console.log(this.chunks);
        console.log('Data saved to disk!') ;
        this._chunks = [] ;
    }
}
