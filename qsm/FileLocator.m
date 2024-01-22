classdef FileLocator
    %FILELOCATOR Locates where files should be
    
    properties
        Dir_top
        Dir_dicoms
    end
    
    methods
        function obj = FileLocator(dir_top)
            %FILELOCATOR Construct an instance of this class
            if ~endsWith(dir_top, '/')
                dir_top = [dir_top '/'];
            end
            obj.Dir_top = dir_top;
            obj.Dir_dicoms = [obj.Dir_top 'dicoms/'];
        end

        function d = GetDicomEchoDir(obj, echoNum)
            d = strcat(obj.Dir_dicoms, 'echo', num2str(echoNum), '/');
        end

        function d = GetDicomDir(obj, echoNum, imgType)
            d = strcat(obj.GetDicomEchoDir(echoNum), imgType, '/');
        end

        function d = GetEchoDir(obj, echoNum)
            d = strcat(obj.Dir_top, 'echo', num2str(echoNum), '/');
        end
        
        function outputArg = GetReal(obj, echoNum)
            %GETREAL Returns path to the real image
            outputArg = obj.GetBasicType(echoNum, 'real');
        end
 
        function outputArg = GetImaginary(obj, echoNum)
            %GetImaginary Returns path to the imaginary image
            outputArg = obj.GetBasicType(echoNum, 'imaginary');
        end

        function outputArg = GetPhase(obj, echoNum)
            %GetPhase Returns path to the phase image
            outputArg = obj.GetBasicType(echoNum, 'phase');
        end

        function outputArg = GetMagnitude(obj, echoNum)
            %GetMagnitude Returns path to the magnitude image
            outputArg = obj.GetBasicType(echoNum, 'mag');
        end

        function outputArg = GetBasicType(obj, echoNum, imgType)
            %GetBasicType Returns path to an image like real or imaginary
            outputArg = obj.GetBasicTypeWithSuffix(echoNum, imgType, "nii");
        end

        function outputArg = GetBasicTypeWithSuffix(obj, echoNum, imgType, suffix)
            %GetBasicType Returns path to an image like real or imaginary
            outputArg = strcat(obj.GetEchoDir(echoNum), imgType, '.', suffix);
        end

        function outputArg = GetMagnitude_AllEchos(obj)
            %GetMagnitude_AllEchos Returns path to the magnitude image for all
            %echoes
            outputArg = strcat(obj.Dir_top, 'mag_allEchos.nii');
        end

        function outputArg = GetPhase_AllEchos(obj)
            %GetPhase_AllEchos Returns path to the phase image for all
            %echoes
            outputArg = strcat(obj.Dir_top, 'phase_allEchos.nii');
        end
    end
end

